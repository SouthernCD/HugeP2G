from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from hugep2g.src.genblasta import genblasta_run, read_genblasta_results
from hugep2g.src.genewise import genewise_run, genewise_output_parser, gff_rec_obj_rename, block_dict_2_gff_rec
from hugep2g.config import wu_blast_path
from interlap import InterLap
from toolbiox.api.common.genome.blast import outfmt6_read_big
from toolbiox.lib.common.genome.genome_feature2 import fancy_name_parse, add_gfs_into_db, read_gff_file, write_gff_file, get_gf_from_db
from toolbiox.lib.common.genome.seq_base import read_fasta_by_faidx, reverse_complement, read_fasta
from toolbiox.lib.common.math.bin_and_window import split_sequence_to_bins
from toolbiox.lib.common.os import have_file, mkdir, rmdir, copy_file, cmd_run, multiprocess_running
from toolbiox.lib.common.util import logging_init
from toolbiox.lib.common.fileIO import tsv_file_dict_parse
from toolbiox.lib.xuyuxing.math.set_operating import overlap_with_interlap_set
import os
# import pandas as pd
import re
import toolbiox.lib.common.sqlite_command as sc

# from toolbiox.api.xuyuxing.genome.repeatmasker import repeatmasker_parser
# from toolbiox.lib.common.math.interval import merge_intervals
#
# def get_skip_file(gff_file, repeat_file, output_file):
#     gene_dict = read_gff_file(gff_file)['gene']
#     repeat_dict = repeatmasker_parser(repeat_file)

#     skiped_range = {}
#     for gene_id in gene_dict:
#         gene = gene_dict[gene_id]

#         if gene.chr_id not in skiped_range:
#             skiped_range[gene.chr_id] = []

#         for mRNA in gene_dict[gene_id].sub_features:
#             for cds in mRNA.sub_features:
#                 if cds.type == 'CDS':
#                     skiped_range[gene.chr_id].append(cds.range)

#     for repeat_type in repeat_dict:
#         for repeat_family_id in repeat_dict[repeat_type]:
#             repeat_family = repeat_dict[repeat_type][repeat_family_id]
#             for repeat_case in repeat_family.case_list:

#                 repeat_case.q_name

#                 if repeat_case.q_name not in skiped_range:
#                     skiped_range[repeat_case.q_name] = []

#                 skiped_range[repeat_case.q_name].append(
#                     (repeat_case.q_start, repeat_case.q_end))

#     for contig in skiped_range:
#         skiped_range[contig] = merge_intervals(skiped_range[contig])

#     with open(output_file, 'w') as f:
#         for contig in skiped_range:
#             for range_tmp in skiped_range[contig]:
#                 f.write("%s\t%d\t%d\n" % (contig, range_tmp[0], range_tmp[1]))


def read_skip_file(skip_file):
    skip_range_dict = {}

    if skip_file:

        with open(skip_file, 'r') as f:
            for each_line in f:

                each_line = each_line.strip()
                chr_id, start, end = each_line.split()
                start = int(start)
                end = int(end)

                if chr_id not in skip_range_dict:
                    skip_range_dict[chr_id] = []

                skip_range_dict[chr_id].append(
                    (min(start, end), max(start, end)))

        for i in skip_range_dict:
            skip_range_dict[i] = InterLap(skip_range_dict[i])

    return skip_range_dict


def init_query_meta_db(query_seq_list, query_fasta_dict, query_sqlite_db, split):
    # remove old db
    # print(query_sqlite_db)
    if os.path.exists(query_sqlite_db):
        cmd_run("rm %s" % query_sqlite_db)

    # build query_seq_table
    work_dir_hash = {}
    dir_name_hash = {}
    for index, start, end in split_sequence_to_bins(len(query_seq_list), split, start=1):
        dir_name = str(start) + "_" + str(end)
        dir_name_hash[dir_name] = []
        query_seq_in_range = query_seq_list[start - 1: end]

        for q_id in query_seq_in_range:
            dir_name_hash[dir_name].append(q_id)
            work_dir_hash[q_id] = dir_name

    sc.init_sql_db(query_sqlite_db, 'query_seq_table', [
                   'id', 'q_id', 'q_len', 'sub_dir'], remove_old_db=True)

    record_list = [(i, query_seq_list[i], len(query_fasta_dict[query_seq_list[i]].faidx), work_dir_hash[query_seq_list[i]])
                   for i in range(len(query_seq_list))]

    sc.sqlite_write(record_list, query_sqlite_db,
                    'query_seq_table', ['id', 'q_id', 'q_len', 'sub_dir'])

    # gba tables
    gba_col = ['id', 'q_id', 'q_start', 'q_end', 's_contig', 's_start', 's_end', 'strand',
               'rank', 'h_name', 'start_e', 'end_e', 'gws_s_base', 'cover', 'skip_cover', 'gws_check']

    table_columns_dict = {}
    for dir_name in dir_name_hash:
        table_columns_dict['work'+dir_name] = gba_col

    sc.init_sql_db_many_table(
        query_sqlite_db, table_columns_dict, remove_old_db=False)


def build_genblasta_dir_tree(query_dir, query_gene_fasta, query_seq_list, split):
    seqdict, seqname_list = read_fasta(query_gene_fasta)
    query_fasta_dict = {i: seqdict[i].seq for i in seqdict}

    for index, start, end in split_sequence_to_bins(len(query_seq_list), split, start=1):
        dir_name = str(start) + "_" + str(end)

        mkdir(query_dir + "/" + dir_name, True)
        query_seq_in_range = query_seq_list[start - 1: end]

        for q_id in query_seq_in_range:
            q_file = query_dir + "/" + dir_name + "/" + q_id + '.fa'

            if os.path.exists(q_file) and os.path.getsize(q_file) != 0:
                continue
            else:
                with open(q_file, 'w') as f:
                    f.write(">%s\n%s\n" % (q_id, query_fasta_dict[q_id]))


def genblasta_results_num(genblasta_results_file):
    cmd_string = "grep \"rank:\" %s | sed \'s/.*rank://g\'" % genblasta_results_file
    flag, output, error = cmd_run(cmd_string, silence=True)

    output_list = [int(i) for i in output.split("\n") if i != '']

    if len(output_list) > 0:
        return max(output_list)
    else:
        return 0


def genblasta_check(gba_file):
    """return: gba_ok, gba_hit_num"""
    if not os.path.exists(gba_file) or os.path.getsize(gba_file) == 0:
        return False, None
    else:
        hit_num = genblasta_results_num(gba_file)
        return True, hit_num


def genblasta_range_not_in_skip_range(gba_file, skip_range_dict, skip_coverage):
    new_range = []

    for q_id, hit_list in read_genblasta_results(gba_file):
        for hit in hit_list:
            query_range = hit.rangeA
            subject_range = hit.rangeB

            skip_flag = False
            if subject_range.chr_id in skip_range_dict:
                cover_length = 0
                total_length = 0

                for hsp in hit.hsp_list:
                    overlap_ratio, overlap_length, overlap_list = overlap_with_interlap_set(
                        hsp.rangeB.range, skip_range_dict[subject_range.chr_id])
                    total_length += hsp.rangeB.len()
                    cover_length += overlap_length

                coverage_ratio = cover_length / total_length

                if coverage_ratio >= skip_coverage:
                    skip_flag = True

            else:
                coverage_ratio = 0.0

            if not skip_flag:
                query_range.get_fancy_name()
                subject_range.get_fancy_name()
                new_range.append(
                    (query_range.fancy_name, subject_range.fancy_name, hit.cover, coverage_ratio))

    return new_range


def get_dir_name_hash(meta_sqlite_db, reverse=False):
    query_gene_info = sc.sqlite_select(
        meta_sqlite_db, 'query_seq_table', ['q_id', 'sub_dir'])
    dir_name_hash = {}

    if reverse:
        return {q_id: sub_dir for q_id, sub_dir in query_gene_info}
    else:
        for q_id, sub_dir in query_gene_info:
            if sub_dir not in dir_name_hash:
                dir_name_hash[sub_dir] = []
            dir_name_hash[sub_dir].append(q_id)

        return dir_name_hash


def get_gws_s_base(contig, start, end, strand, contig_length, expand_range_ratio=0.1):

    expand_range = int(
        (abs(end - start) + 1) * expand_range_ratio)

    start_e = max(
        1, min(start, end) - expand_range)
    end_e = min(contig_length, max(start, end) + expand_range)

    s_base = start_e
    if strand == "-":
        s_base = end_e

    return s_base, start_e, end_e


def genblasta_dir_results_check(query_dir, target_fasta_file, skip_range_file, num_threads, log_file, skip_coverage, genblasta_hit_num):
    annotated_range_dict = {}

    query_sqlite_db = query_dir + "/genblasta.db"

    dir_name_hash = get_dir_name_hash(query_sqlite_db)

    for sub_dir in dir_name_hash:

        # check gba

        parsed_query_gene_list = list(
            set([i[0] for i in sc.sqlite_select(query_sqlite_db, 'work'+sub_dir, ['q_id'])]))
        unparsed_query_gene_list = list(
            set(dir_name_hash[sub_dir]) - set(parsed_query_gene_list))

        check_args_list = []
        for q_id in unparsed_query_gene_list:
            gba_file = query_dir + "/" + sub_dir + "/" + q_id + ".fa.gba"
            check_args_list.append((gba_file,))

        check_out = multiprocess_running(genblasta_check, check_args_list, num_threads,
                                         log_file=log_file, args_id_list=unparsed_query_gene_list)

        check_ok_list = [i for i in check_out if check_out[i]['output'][0]]

        # compare gff

        compare_args_list = []

        if len(check_ok_list) > 0 and len(annotated_range_dict) == 0:
            skip_range_dict = read_skip_file(skip_range_file)

        # for q_id in check_ok_list:
        #     gba_file = query_dir + "/" + sub_dir + "/" + q_id + ".fa.gba"
        #     compare_args_list.append(
        #         (gba_file, skip_range_dict, skip_coverage))

        # compare_out = multiprocess_running(genblasta_range_not_in_skip_range, compare_args_list,
        #                                    num_threads, log_file=log_file, args_id_list=check_ok_list)

        compare_out = {}
        for q_id in check_ok_list:
            gba_file = query_dir + "/" + sub_dir + "/" + q_id + ".fa.gba"
            compare_out[q_id] = {}
            compare_out[q_id]['output'] = genblasta_range_not_in_skip_range(
                gba_file, skip_range_dict, skip_coverage)

        # write sqlite

        record_id_now = sc.sql_table_row_num(query_sqlite_db, 'work'+sub_dir)
        target_fasta_dict = read_fasta_by_faidx(target_fasta_file)

        waitting_for_load = []
        for q_id in check_ok_list:

            if len(compare_out[q_id]['output']) == 0:
                waitting_for_load.append(
                    (record_id_now, q_id, 0, 0, "None", 0, 0, ".", 0, "None", 0, 0, 0, 0.0, 0.0, True))
                record_id_now += 1
            else:
                h_num = 0
                for q_range, h_range, cover, skip_cover in compare_out[q_id]['output']:
                    s_contig, s_start, s_end, s_strand = fancy_name_parse(
                        h_range)
                    q_contig, q_start, q_end, q_strand = fancy_name_parse(
                        q_range)

                    h_name = q_id + "_" + str(h_num)

                    s_base, start_e, end_e = get_gws_s_base(
                        s_contig, s_start, s_end, s_strand, len(target_fasta_dict[s_contig].faidx), 0)

                    waitting_for_load.append(
                        (record_id_now, q_id, q_start, q_end, s_contig, s_start, s_end, s_strand, h_num, h_name, start_e, end_e, s_base, cover, skip_cover, False))
                    h_num += 1
                    record_id_now += 1

                    if h_num >= genblasta_hit_num:
                        break

            if len(waitting_for_load) > 10000:
                sc.sqlite_write(waitting_for_load, query_sqlite_db, 'work'+sub_dir, [
                                'id', 'q_id', 'q_start', 'q_end', 's_contig', 's_start', 's_end', 'strand', 'rank', 'h_name', 'start_e', 'end_e', 'gws_s_base', 'cover', 'skip_cover', 'gws_check'])
                waitting_for_load = []

        if not len(waitting_for_load) == 0:
            sc.sqlite_write(waitting_for_load, query_sqlite_db, 'work'+sub_dir, [
                            'id', 'q_id', 'q_start', 'q_end', 's_contig', 's_start', 's_end', 'strand', 'rank', 'h_name', 'start_e', 'end_e', 'gws_s_base', 'cover', 'skip_cover', 'gws_check'])

        del waitting_for_load


def seq_stats(meta_sqlite_db, check_key=None, return_key='q_id'):
    table_list = list(sc.check_sql_table(meta_sqlite_db))
    table_list.remove('query_seq_table')

    big_output_list = []

    for table_name in table_list:

        if check_key is None:
            output_list = sc.sqlite_select(
                meta_sqlite_db, table_name, column_list=[return_key])
        else:
            output_list = sc.sqlite_select(meta_sqlite_db, table_name, column_list=[
                return_key], key_name=check_key, value_tuple=(True,))

        output_list = list(set([i[0] for i in output_list]))

        big_output_list.extend(output_list)

    return list(set(big_output_list))


def read_gws_file(gws_file, query_length, strand, gws_base, h_name, q_name, s_name):

    if not os.path.exists(gws_file) or os.path.getsize(gws_file) == 0:
        return False, None, None, None

    try:

        gff_dict = {}
        pro_dict = {}
        cds_dict = {}

        block_dict, gff_output, block_pro_seq, block_cds_seq = genewise_output_parser(
            gws_file)

        # stat
        # [re.sub(' ', '', a[0][i]['block'][0]) for i in a[0]]
        block_dict

        aln_len = sum([len(re.sub(' ', '', block_dict[i]['block'][0]))
                       for i in block_dict])

        pst_len = sum([len(re.sub(' ', '', block_dict[i]['block'][1]))
                       for i in block_dict])

        idt_len = sum([len(re.sub('[ +]', '', block_dict[i]['block'][1]))
                       for i in block_dict])

        contig_rec = block_dict_2_gff_rec(block_dict, gff_output, query_length)

        if strand == "+":
            strand = 1
        else:
            strand = -1

        renamed_contig_rec = gff_rec_obj_rename(
            contig_rec, hit_name=h_name, q_name=q_name, s_name=s_name, s_base=gws_base, strand=strand)

        rec = SeqRecord(Seq(""), s_name)
        gff_dict[s_name] = rec

        gff_dict[s_name].features.append(
            renamed_contig_rec.features[0])

        pro_dict[h_name] = block_pro_seq
        cds_dict[h_name] = block_cds_seq

        gff_file = gws_file + ".gff3"

        rec_obj_list = [gff_dict[i] for i in gff_dict]
        with open(gff_file, "w") as f:
            GFF.write(rec_obj_list, f)

        gff_dict = read_gff_file(gff_file)

        for i in gff_dict:
            for j in gff_dict[i]:
                gf = gff_dict[i][j]

        gf.sub_features[0].qualifiers['aln_len'] = aln_len
        gf.sub_features[0].qualifiers['pst_len'] = pst_len
        gf.sub_features[0].qualifiers['idt_len'] = idt_len

        os.remove(gff_file)

        return True, gf, cds_dict, pro_dict

    except:

        return False, None, None, None


def get_q_len_from_db(query_sqlite_db):
    output_dict = {}
    for q_id, q_len in sc.sqlite_select(query_sqlite_db, 'query_seq_table', ['q_id', 'q_len']):
        output_dict[q_id] = q_len
    return output_dict


def genewise_dir_results_check(query_dir, query_sqlite_db, sub_dir, q_len_dict, output_gws_db, output_pep_file, output_cds_file, num_threads, log_file):

    args_list = []
    id_list = []

    for (id_tmp, q_id, q_start, q_end, s_contig, start, end, strand, rank, h_name, start_e, end_e, gws_s_base, cover, skip_cover, gws_check) in sc.sqlite_select(query_sqlite_db, 'work'+sub_dir, None, 'gws_check', (0,)):
        output_gws_file = query_dir + "/" + sub_dir + \
            "/" + q_id + "/" + str(rank) + "/output.gws"
        if os.path.exists(output_gws_file) and os.path.getsize(output_gws_file) > 0:
            id_list.append(id_tmp)
            args_list.append((output_gws_file, q_len_dict[q_id], strand, gws_s_base,
                              h_name, q_id, s_contig))

    if len(args_list) > 0:
        check_out = multiprocess_running(
            read_gws_file, args_list, num_threads, log_file=log_file, args_id_list=id_list)

        gf_big_list = []
        for id_tmp in check_out:
            if check_out[id_tmp]['output'][0]:
                flag, gf, cds_dict, pro_dict = check_out[id_tmp]['output']
                if flag:
                    gf_big_list.append(gf)
                    with open(output_cds_file, 'a') as f:
                        for i in cds_dict:
                            f.write(">%s\n%s\n" % (i, cds_dict[i]))

                    with open(output_pep_file, 'a') as f:
                        for i in pro_dict:
                            f.write(">%s\n%s\n" % (i, pro_dict[i]))

        if len(gf_big_list) > 0:
            add_gfs_into_db(gf_big_list, output_gws_db)

            gws_check_out_list = [(id_tmp, check_out[id_tmp]['output'][0])
                                  for id_tmp in check_out if check_out[id_tmp]['output'][0]]

            sc.sqlite_update(query_sqlite_db, 'work'+sub_dir,
                             'id', 'gws_check', gws_check_out_list)


def prepare_genewise_running(query_sqlite_db, query_dir, sub_dir, query_seq_dict, target_genome_dict):

    genewise_args_list = []

    for (id_tmp, q_id, q_start, q_end, s_contig, start, end, strand, rank, h_name, start_e, end_e, gws_s_base, cover, skip_cover, gws_check) in sc.sqlite_select(query_sqlite_db, 'work'+sub_dir, None, 'gws_check', (0,)):
        query_gws_dir = query_dir + "/" + sub_dir + "/" + q_id + "/"
        mkdir(query_gws_dir, True)
        query_rank_gws_dir = query_gws_dir + "/" + str(rank) + "/"
        mkdir(query_rank_gws_dir, True)

        seq_str = target_genome_dict[s_contig][start_e - 1: end_e].upper()

        if strand == "-":
            seq_str = reverse_complement(seq_str)

        if not (os.path.exists(query_rank_gws_dir + "/subject.fa") and os.path.getsize(
                query_rank_gws_dir + "/subject.fa") != 0):
            with open(query_rank_gws_dir + "/subject.fa", 'w') as f:
                f.write(">subject\n%s\n" % seq_str)
        if not (os.path.exists(query_rank_gws_dir + "/query.fa") and os.path.getsize(
                query_rank_gws_dir + "/query.fa") != 0):
            with open(query_rank_gws_dir + "/query.fa", 'w') as f:
                f.write(">query\n%s\n" % query_seq_dict[q_id])

        genewise_args_list.append((query_rank_gws_dir + "/query.fa", query_rank_gws_dir +
                                   "/subject.fa", query_rank_gws_dir + "/output.gws", '+', True, q_start, q_end))

    return genewise_args_list


def parse_diamond_out(diamond_file):

    ev_diamond_info = {}
    for query_record in outfmt6_read_big(diamond_file):
        q_id = query_record.qID
        ev_diamond_info[q_id] = {}

        q_parent_id = re.match('^(.*)_\d+$', q_id).groups()[0]

        for hit_record in query_record.hit:
            hit_id = hit_record.Hit_id
            if hit_id == q_parent_id:
                ev_diamond_info[q_id]["hit_rank"] = hit_record.Hit_num
                top_hsp = hit_record.hsp[0]
                ev_diamond_info[q_id]["identity"] = top_hsp.Hsp_identical_ratio
                ev_diamond_info[q_id]["bit_score"] = top_hsp.Hsp_bit_score
                ev_diamond_info[q_id]["evalue"] = top_hsp.Hsp_evalue
                ev_diamond_info[q_id]["align_len"] = top_hsp.Hsp_align_len

    return ev_diamond_info


def merge_all_info(query_dir):
    query_fasta = os.path.join(query_dir, 'query.protein.fasta')
    # genblasta_db = os.path.join(query_dir, 'genblasta.db')
    genewise_db = os.path.join(query_dir, 'genewise.db')
    diamond_db = os.path.join(query_dir, 'genewise.pep.bls.db')
    pep_fasta = os.path.join(query_dir, 'genewise.pep.fa')
    cds_fasta = os.path.join(query_dir, 'genewise.cds.fa')

    gf_dict = get_gf_from_db(genewise_db)
    diamond_dict = sc.retrieve_dict_from_db(diamond_db)

    pep_seq_dict = read_fasta_by_faidx(pep_fasta)
    cds_seq_dict = read_fasta_by_faidx(cds_fasta)

    query_seq_dict = read_fasta_by_faidx(query_fasta)
    for i in query_seq_dict:
        query_seq_dict[i].seq = str(query_seq_dict[i].seq)
        query_seq_dict[i].seq_type = 'prot'
        query_seq_dict[i].subjects = []
        query_seq_dict[i].low_quality_subject = []
        del query_seq_dict[i].faidx

    for gf_id in gf_dict:
        gf = gf_dict[gf_id]
        mgf = gf.sub_features[0]
        mgf.qualifiers['AA'] = str(pep_seq_dict[gf_id].seq)
        mgf.qualifiers['CDS'] = str(cds_seq_dict[gf_id].seq)
        query_seq_id = gf.sub_features[0].qualifiers['Target'][0]

        if gf_id in diamond_dict and 'hit_rank' in diamond_dict[gf_id]:
            # mgf.qualifiers['bls_rank'] = int(diamond_dict[gf_id]['hit_rank'])
            mgf.qualifiers['evalue'] = float(diamond_dict[gf_id]['evalue'])
            query_seq_dict[query_seq_id].subjects.append(gf)
        else:
            # mgf.qualifiers['bls_rank'] = -1
            mgf.qualifiers['evalue'] = 999
            query_seq_dict[query_seq_id].low_quality_subject.append(gf)

    for query_seq_id in query_seq_dict:
        subjects_list = query_seq_dict[query_seq_id].subjects
        subjects_list = sorted(
            subjects_list, key=lambda x: x.sub_features[0].qualifiers['evalue'], reverse=False)
        query_seq_dict[query_seq_id].subjects = subjects_list

    return query_seq_dict


def build_hugep2g(args):

    # init project
    # print arguments to log file
    args.work_dir = os.path.abspath(args.work_dir)
    mkdir(args.work_dir, True)

    args.log_file = args.work_dir + "/log"

    logger = logging_init("Huge Proteins Map To Genome", args.log_file)

    args_info_string = "Argument list:\n"
    attrs = vars(args)

    for item in attrs.items():
        args_info_string = args_info_string + ("%s: %s\n" % item)

    logger.info(args_info_string)

    # read query info table

    info_dict = tsv_file_dict_parse(args.query_protein_table)
    query_file_dict = {}
    for i in info_dict:
        query_file_dict[info_dict[i]['sp_id']] = info_dict[i]['pt_file']

    # query_db = pd.read_excel(args.query_protein_table)
    # query_file_dict = {}
    # for index in query_db.index:
    #     query_id = str(query_db.loc[index]['sp_id'])
    #     query_pt_file = str(query_db.loc[index]['pt_file'])

    #     if not pd.isna(query_id) and not pd.isna(query_pt_file):
    #         query_file_dict[query_id] = query_pt_file

    # diamond_dir = args.work_dir + "/diamond_dir"
    # huge_ref_pt_file = diamond_dir + "/huge_ref_pt.fasta"
    # diamond_db = huge_ref_pt_file + ".dmnd"
    # if args.force_redo or not os.path.exists(diamond_dir):
    #     mkdir(diamond_dir, True)

    #     for query_id in query_file_dict:
    #         pt_file = query_file_dict[query_id]
    #         cmd_string = "cat %s >> %s" % (pt_file, huge_ref_pt_file)
    #         cmd_run(cmd_string, silence=True)

    #     cmd_string = "diamond makedb --in %s --db %s" % (
    #         huge_ref_pt_file, huge_ref_pt_file)
    #     cmd_run(cmd_string, silence=True)

    printer_string = ""
    for i in query_file_dict:
        printer_string += "ID: %s\tfile: %s\n" % (i, query_file_dict[i])
    logger.info("There are %d query genome files: \n" %
                len(query_file_dict) + printer_string)

    # copy target_genome_fasta to work_dir
    logger.info("copy target genome to work dir")
    args.target_genome_fasta = copy_file(
        args.target_genome_fasta, args.work_dir, True)

    # copy target_genome_gff to work_dir
    logger.info("copy skip_range_file to work dir\n")
    if args.skip_range_file:
        args.skip_range_file = copy_file(
            args.skip_range_file, args.work_dir, True)

    # build dir tree
    logger.info("Build work dir tree\n")

    run_dir = args.work_dir + "/work_dir"

    mkdir(run_dir, True)

    for i in query_file_dict:
        mkdir(run_dir + "/" + i, True)
        query_file_dict[i] = copy_file(
            query_file_dict[i], os.path.join(run_dir, i, 'query.protein.fasta'), True)

        diamond_db = query_file_dict[i] + ".dmnd"
        if args.force_redo or not have_file(diamond_db):

            cmd_string = "diamond makedb --in %s --db %s" % (
                query_file_dict[i], query_file_dict[i])
            cmd_run(cmd_string, silence=True)

    # genblasta
    logger.info("Genblasta: ")

    for query_speci_id in query_file_dict:
        logger.info("prepare on %s" % query_speci_id)

        query_fasta_dict = read_fasta_by_faidx(query_file_dict[query_speci_id])
        query_seq_list = list(query_fasta_dict.keys())

        query_dir = "%s/%s" % (run_dir, query_speci_id)

        query_sqlite_db = "%s/%s/genblasta.db" % (run_dir, query_speci_id)

        # create sqlite db
        if args.force_redo or not os.path.exists(query_sqlite_db):
            init_query_meta_db(
                query_seq_list, query_fasta_dict, query_sqlite_db, args.seq_num_in_subdir)

        # build dir tree split sequence to file
        dir_file_list = os.listdir(query_dir)
        ok_flag = True
        for i in [str(i[1])+"_"+str(i[2]) for i in split_sequence_to_bins(len(query_seq_list), args.seq_num_in_subdir, start=1)]:
            if i not in dir_file_list:
                ok_flag = False

        if args.force_redo or ok_flag is False:
            build_genblasta_dir_tree(
                query_dir, query_file_dict[query_speci_id], query_seq_list, args.seq_num_in_subdir)

    # running genblasta
    logger.info("running genblasta\n")
    cmd_string = "%s/xdformat -n %s" % (wu_blast_path,
                                        args.target_genome_fasta)
    tmp_out = cmd_run(cmd_string, silence=True)

    for query_speci_id in query_file_dict:
        logger.info("\tparse %s" % query_speci_id)

        query_dir = "%s/%s" % (run_dir, query_speci_id)
        meta_sqlite_db = query_dir+"/genblasta.db"

        query_fasta_dict = read_fasta_by_faidx(query_file_dict[query_speci_id])
        query_seq_list = list(query_fasta_dict.keys())

        genblasta_dir_results_check(query_dir, args.target_genome_fasta,
                                    args.skip_range_file, args.num_threads, None, args.skip_coverage, args.genblasta_hit_num)

        gba_passed_id_list = seq_stats(meta_sqlite_db)

        un_gba_query_gene_list = list(
            set(query_seq_list) - set(gba_passed_id_list))

        logger.info("\tThere are %s genes: %s already have good hits, %s need run now" % (
            len(query_seq_list), len(gba_passed_id_list), len(un_gba_query_gene_list)))

        dir_name_hash = get_dir_name_hash(meta_sqlite_db, True)

        args_list = []
        for q_id in un_gba_query_gene_list:
            q_file = query_dir + "/" + dir_name_hash[q_id] + "/" + q_id + ".fa"
            gba_file = q_file + ".gba"
            args_list.append((q_file, args.target_genome_fasta,
                              gba_file, args.gene_coverage, args.genblasta_hit_num))

        genblasta_out = multiprocess_running(
            genblasta_run, args_list, args.num_threads, log_file=args.log_file, timeout=1800)

        genblasta_dir_results_check(query_dir, args.target_genome_fasta,
                                    args.skip_range_file, args.num_threads, None, args.skip_coverage, args.genblasta_hit_num)

        gba_passed_id_list = seq_stats(meta_sqlite_db)
        un_gba_query_gene_list = list(
            set(query_seq_list) - set(gba_passed_id_list))

        logger.info("After genblasta: %s genes have good hits, %s failed" % (
            len(gba_passed_id_list), len(un_gba_query_gene_list)))

    logger.info("\nGeneWise: ")

    target_genome_dict = read_fasta_by_faidx(args.target_genome_fasta)
    target_genome_dict = {
        i: str(target_genome_dict[i].seq) for i in target_genome_dict}

    for query_speci_id in query_file_dict:
        query_dir = "%s/%s" % (run_dir, query_speci_id)
        logger.info("\nparse %s" % query_speci_id)
        query_seq_dict = read_fasta_by_faidx(query_file_dict[query_speci_id])
        query_seq_dict = {i: str(query_seq_dict[i].seq)
                          for i in query_seq_dict}

        query_sqlite_db = query_dir + "/genblasta.db"
        dir_name_hash = get_dir_name_hash(query_sqlite_db)
        dir_name_hash_reverse = get_dir_name_hash(query_sqlite_db, True)
        q_len_dict = get_q_len_from_db(query_sqlite_db)

        output_gws_db = query_dir + "/genewise.db"
        output_pep_file = query_dir + "/genewise.pep.fa"
        output_cds_file = query_dir + "/genewise.cds.fa"

        for sub_dir in dir_name_hash:
            logger.info("\tparse sub_dir: %s" % sub_dir)

            # # parse old genewise output
            genewise_dir_results_check(query_dir, query_sqlite_db, sub_dir, q_len_dict, output_gws_db,
                                       output_pep_file, output_cds_file, args.num_threads, None)

            # prepare genewise input file
            genewise_args_list = prepare_genewise_running(
                query_sqlite_db, query_dir, sub_dir, query_seq_dict, target_genome_dict)

            # running genewise
            multiprocess_running(
                genewise_run, genewise_args_list, args.num_threads, log_file=args.log_file, timeout=1800)

            # reparse genewise output
            logger.info("\tread results in sub_dir: %s" % sub_dir)
            genewise_dir_results_check(query_dir, query_sqlite_db, sub_dir, q_len_dict, output_gws_db,
                                       output_pep_file, output_cds_file, args.num_threads, None)

            rmdir(os.path.join(query_dir, sub_dir))

        # make diamond for pep seq
        logger.info("running diamond")

        diamond_db = query_file_dict[query_speci_id] + ".dmnd"
        diamond_out = query_dir + "/genewise.pep.bls"
        diamond_out_db = query_dir + "/genewise.pep.bls.db"

        if not have_file(diamond_out_db):
            cmd_string = "diamond blastp -q %s -k 100 -e 1e-5 -d %s -o %s -f 6 -p %d" % (
                output_pep_file, diamond_db, diamond_out, args.num_threads)
            cmd_run(cmd_string, silence=True)

            ev_diamond_info = parse_diamond_out(diamond_out)
            sc.store_dict_to_db(ev_diamond_info, diamond_out_db)

        # merge_info
        logger.info("make output")
        out_db = query_dir + "/output.db"
        merged_info_dict = merge_all_info(query_dir)
        sc.store_dict_to_db(merged_info_dict, out_db)

        gf_list = []
        for q_id in merged_info_dict:
            q_record = merged_info_dict[q_id]
            gf_list.extend(q_record.subjects + q_record.low_quality_subject)

        for gf in gf_list:
            mgf = gf.sub_features[0]
            del mgf.qualifiers['AA']
            del mgf.qualifiers['CDS']

        write_gff_file(gf_list, query_dir + "/genewise.gff3")

        try:
            rmdir(query_dir + "/query.protein.fasta")
            rmdir(query_dir + "/query.protein.fasta.fai")
            rmdir(query_dir + "/query.protein.fasta.dmnd")
            rmdir(query_dir + "/genewise.db")
            rmdir(query_dir + "/genblasta.db")
            rmdir(query_dir + "/genewise.pep.bls")
            rmdir(query_dir + "/genewise.pep.bls.db")
            rmdir(query_dir + "/genewise.pep.fa.fai")
            rmdir(query_dir + "/genewise.cds.fa.fai")
        except:
            pass

        logger.info("Finished %s\n" % query_speci_id)

# ##################
# def extract_one_evidences(id_tmp, gf_name, type_tmp, contig_name, start, end, strand, daughter, qualifiers, gws_db, OG_support="", taxon_id="", table_name=""):
#     gf = gf_info_retrieval((id_tmp, gf_name, type_tmp, contig_name,
#                             start, end, strand, daughter, qualifiers), 0, 'A', gws_db)

#     mRNA_gf = gf.sub_features[0]
#     qualifiers = mRNA_gf.qualifiers

#     # query coverage
#     coverage = min((int(qualifiers['Target_End'][0]) - int(
#         qualifiers['Target_Start'][0]) + 1) / int(qualifiers['Target_Length'][0]), 1.00)

#     # aln aa
#     aln_aa_len = int(
#         qualifiers['Target_End'][0]) - int(qualifiers['Target_Start'][0]) + 1

#     # identity
#     identity = int(qualifiers['idt_len']) / \
#         int(qualifiers['aln_len'])

#     q_id = re.match('^(.*)_\d+$', gf_name).groups()[0]

#     gf.evidence_indicator = {
#         "taxon_id": taxon_id,
#         "query_id": q_id,
#         "query_coverage": coverage,
#         "aln_aa_len": aln_aa_len,
#         "identity": identity,
#         "query_length": int(qualifiers['Target_Length'][0]),
#         "score": max(0.1, float(gf.sub_features[0].qualifiers['score'][0])),
#         "OG_support": OG_support
#     }

#     gf.score_dict = {
#         "score": max(0.1, float(gf.sub_features[0].qualifiers['score'][0]))
#     }

#     gf.db_path = "%s.%s.%s" % (taxon_id, table_name, id_tmp)

#     return gf


# def extract_gf_id(WPGmapper_dir):
#     WPGmapper_dict = {'genblasta': os.path.join(WPGmapper_dir, "genblasta.db"),
#                       'genewise': os.path.join(WPGmapper_dir, "genewise.db"),
#                       'pep_db': os.path.join(WPGmapper_dir, "genewise.pep.db"),
#                       'cds_db': os.path.join(WPGmapper_dir, "genewise.cds.db"),
#                       'diamond': os.path.join(WPGmapper_dir, "genewise.pep.bls")
#                       }

#     gws_db = WPGmapper_dict['genewise']
#     meta_dict = get_gf_db_meta_dict(gws_db)

#     gf_id_list = []
#     for table_name in meta_dict['gene'][0]:
#         # print(table_name)
#         sql_cmd_string = 'SELECT * FROM %s' % table_name

#         data_list = sc.sqlite_execute(sql_cmd_string, gws_db)

#         for gf_db_info_tuple in data_list:
#             id_tmp, gf_name, type_tmp, contig_name, start, end, strand, daughter, qualifiers = gf_db_info_tuple

#             q_id = gf_name.split("_")[0]

#             gf_id_list.append(gf_name)

#     return gf_id_list


# def extract_all_evidences(ref_WPGmapper_dict, evidence_id_list=None, seq_support_dict=None, num_thread=20, log_file=None):
#     module_log = logging_init("extract_all_evidences", log_file)

#     # get all useable evidence

#     module_log.info('\tBegin: extract parent diamond blast results')
#     diamond_args_list = []

#     for taxon_id in ref_WPGmapper_dict:
#         diamond_file = ref_WPGmapper_dict[taxon_id]['diamond']
#         diamond_args_list.append((diamond_file, ))

#     diamond_output = multiprocess_running(
#         parse_diamond_out, diamond_args_list, num_thread, silence=True, log_file=log_file)

#     diamond_seq_info = {}
#     num = 0
#     for i in diamond_output:
#         for j in diamond_output[i]['output']:
#             num += 1
#             if len(diamond_output[i]['output'][j]) > 0:
#                 diamond_seq_info[j] = diamond_output[i]['output'][j]

#     module_log.info('\tEnd: extract parent diamond blast results')

#     module_log.info('\tBegin: extract evidence gf from database')

#     if not evidence_id_list is None:
#         evidence_id_dict = {i: 0 for i in evidence_id_list}

#     ev_args_list = []
#     ev_args_id_list = []

#     for taxon_id in ref_WPGmapper_dict:
#         # print(taxon_id)

#         gws_db = ref_WPGmapper_dict[taxon_id]['genewise']
#         meta_dict = get_gf_db_meta_dict(gws_db)

#         # read mRNA info

#         for table_name in meta_dict['gene'][0]:
#             # print(table_name)
#             sql_cmd_string = 'SELECT * FROM %s' % table_name

#             data_list = sc.sqlite_execute(sql_cmd_string, gws_db)
#             for gf_db_info_tuple in data_list:
#                 id_tmp, gf_name, type_tmp, contig_name, start, end, strand, daughter, qualifiers = gf_db_info_tuple

#                 if not evidence_id_list is None:
#                     if gf_name not in evidence_id_dict:
#                         continue

#                 # if gf_name not in diamond_seq_info:
#                 #     continue

#                 q_id = gf_name.split("_")[0]

#                 OG_support = (0, 0)
#                 if not seq_support_dict is None:
#                     if q_id in seq_support_dict:
#                         OG_support = seq_support_dict[q_id]

#                 ev_args_list.append((id_tmp, gf_name, type_tmp, contig_name, start, end,
#                                      strand, daughter, qualifiers, gws_db, OG_support, taxon_id, table_name))
#                 ev_args_id_list.append(gf_name)

#     ev_output = multiprocess_running(
#         extract_one_evidences, ev_args_list, num_thread, silence=True, args_id_list=ev_args_id_list, log_file=log_file)

#     evidence_dict_tmp = OrderedDict()

#     for i in ev_output:
#         gf = ev_output[i]['output']
#         gf.parent_blast = {"hit_rank": 0, "identity": 0,
#                            "bit_score": 0, "evalue": 100}
#         gf.score_dict["bit_score"] = 0.0
#         if i in diamond_seq_info:
#             gf.parent_blast = diamond_seq_info[i]
#             gf.score_dict["bit_score"] = diamond_seq_info[i]["bit_score"]
#         evidence_dict_tmp[i] = gf

#     if not evidence_id_list is None:
#         evidence_dict = OrderedDict()
#         for i in evidence_id_list:
#             if i in evidence_dict_tmp:
#                 evidence_dict[i] = evidence_dict_tmp[i]
#     else:
#         evidence_dict = evidence_dict_tmp

#     module_log.info('\tEnd: extract evidence gf from database')

#     # get aa cds seq
#     module_log.info('\tBegin: extract AA and cds sequence for each gf')
#     species_tmp_dict = {}
#     for i in evidence_dict:
#         taxon_id = evidence_dict[i].evidence_indicator['taxon_id']
#         if taxon_id not in species_tmp_dict:
#             species_tmp_dict[taxon_id] = []
#         species_tmp_dict[taxon_id].append(i)

#     gf_aa_dict = {}
#     gf_cds_dict = {}
#     for taxon_id in species_tmp_dict:
#         gf_aa_dict[taxon_id] = extract_seq_from_sqlite(
#             species_tmp_dict[taxon_id], ref_WPGmapper_dict[taxon_id]['pep_db'])
#         gf_cds_dict[taxon_id] = extract_seq_from_sqlite(
#             species_tmp_dict[taxon_id], ref_WPGmapper_dict[taxon_id]['cds_db'])

#     for i in evidence_dict:
#         gf = evidence_dict[i]
#         taxon_id = gf.evidence_indicator['taxon_id']
#         gf.model_aa_seq = gf_aa_dict[taxon_id][i]
#         gf.model_cds_seq = gf_cds_dict[taxon_id][i]
#         evidence_dict[i] = gf

#     del module_log.handlers[:]

#     return evidence_dict


# def load_map_files(ref_taxon_id_list, WPGmapper_dir):

#     # load map database
#     ref_WPGmapper_dict = {}
#     for taxon_id in ref_taxon_id_list:
#         gba_db_file = WPGmapper_dir + "/WPGmapper/" + taxon_id + "/genblasta.db"
#         gw_db_file = WPGmapper_dir + "/WPGmapper/" + taxon_id + "/genewise.db"
#         gw_pep_file = WPGmapper_dir + "/WPGmapper/" + taxon_id + "/genewise.pep.db"
#         gw_cds_file = WPGmapper_dir + "/WPGmapper/" + taxon_id + "/genewise.cds.db"
#         gw_pep_bls = WPGmapper_dir + "/WPGmapper/" + taxon_id + "/genewise.pep.bls"

#         if have_file(gba_db_file) and have_file(gw_db_file) and have_file(gw_pep_file) and have_file(gw_cds_file) and have_file(gw_pep_bls):
#             ref_WPGmapper_dict[taxon_id] = {}
#             ref_WPGmapper_dict[taxon_id]['genblasta'] = gba_db_file
#             ref_WPGmapper_dict[taxon_id]['genewise'] = gw_db_file
#             ref_WPGmapper_dict[taxon_id]['pep_db'] = gw_pep_file
#             ref_WPGmapper_dict[taxon_id]['cds_db'] = gw_cds_file
#             ref_WPGmapper_dict[taxon_id]['diamond'] = gw_pep_bls

#     return ref_WPGmapper_dict


# def extract_all_evidences_from_one_gws_db(genewise_db):
#     ev_args_list = []
#     ev_args_id_list = []

#     meta_dict = get_gf_db_meta_dict(genewise_db)

#     # read mRNA info

#     for table_name in meta_dict['gene'][0]:
#         # print(table_name)
#         sql_cmd_string = 'SELECT * FROM %s' % table_name

#         data_list = sc.sqlite_execute(sql_cmd_string, genewise_db)
#         for gf_db_info_tuple in data_list:
#             id_tmp, gf_name, type_tmp, contig_name, start, end, strand, daughter, qualifiers = gf_db_info_tuple

#             q_id = gf_name.split("_")[0]

#             ev_args_list.append((id_tmp, gf_name, type_tmp, contig_name, start, end,
#                                  strand, daughter, qualifiers, genewise_db))
#             ev_args_id_list.append(gf_name)

#     ev_output = multiprocess_running(
#         extract_one_evidences, ev_args_list, 56, silence=False, args_id_list=ev_args_id_list)

#     output_dict = {}
#     for i in ev_output:
#         output_dict[i] = ev_output[i]['output']

#     return output_dict


# def get_rep_gff(genewise_db, gff_file):
#     gf_dict = extract_all_evidences_from_one_gws_db(genewise_db)

#     query_gf_dict = {}
#     for gf_id in gf_dict:
#         gf = gf_dict[gf_id]
#         q_id = gf.evidence_indicator['query_id']
#         if not q_id in query_gf_dict:
#             query_gf_dict[q_id] = []
#         query_gf_dict[q_id].append(gf)

#     best_gf_list = []
#     for q_id in query_gf_dict:
#         gf_list = query_gf_dict[q_id]
#         best_gf = max(gf_list, key=lambda x: x.evidence_indicator['score'])
#         best_gf.qualifiers['query'] = q_id
#         best_gf.qualifiers['q_cov'] = best_gf.evidence_indicator['query_coverage'] * 100
#         best_gf.qualifiers['q_idt'] = best_gf.evidence_indicator['identity'] * 100

#         if 'Note' in best_gf.sub_features[0].qualifiers:
#             best_gf.qualifiers['Note'] = 'pseudogene'

#         best_gf_list.append(best_gf)
#     write_gff_file(best_gf_list, gff_file)
if __name__ == "__main__":
    pass