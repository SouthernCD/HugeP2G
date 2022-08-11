from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from hugep2g.config import tmp_work_dir, genblasta_path, genblasta_wu_script
from hugep2g.src.genewise import genewise_run, genewise_output_parser, gff_rec_obj_rename, block_dict_2_gff_rec
from pyfaidx import Fasta
from toolbiox.lib.common.genome.genome_feature2 import HomoPredictResults, Gene, ChrLoci, SimpleRangePair, BlastHspRecord
from toolbiox.lib.common.genome.seq_base import reverse_complement, read_fasta, read_fasta_big, read_fasta_by_faidx
from toolbiox.lib.common.os import mkdir, rmdir, cmd_run
import os
import re
import shutil
import uuid


# def read_genblasta_results(gen_file):
#     with open(gen_file, 'r') as f:
#         query_record = {}
#         query_record['query'] = ""
#         query_record['subject'] = {}
#         for each_line in f:
#             each_line = each_line.strip()
#             # record start
#             matchObj = re.match("^\/\/.*START.*\/\/$", each_line)
#             if matchObj:
#                 continue
#
#             # record query name
#             matchObj = re.match("\/\/for query:\s+(\S+)\s+\/\/", each_line)
#             if matchObj:
#                 query_record['query'] = matchObj.groups()[0]
#                 continue
#
#             # subject record
#             matchObj = re.match(
#                 "^(.*)\|(.*):(\d+)\.\.(\d+)\|(\S+)\|gene cover:(\d+)\((.*)%\)\|score:(.*)\|rank:(\d+)$",
#                 each_line)
#             if matchObj:
#                 q_name, s_name, s_start, s_end, strand, aln, cover, score, rank = matchObj.groups()
#                 query_record['subject'][len(query_record['subject'])] = (
#                     q_name, s_name, s_start, s_end, strand, aln, cover, score, rank)
#                 continue
#
#             # record end
#             matchObj = re.match("^\/\/.*END.*\/\/$", each_line)
#             if matchObj:
#                 yield query_record
#                 query_record = {}
#                 query_record['query'] = ""
#                 query_record['subject'] = {}
#                 continue

class GenBlastARecord(SimpleRangePair):
    def __init__(self, rangeA, rangeB, score, aln, cover, rank, hsp_list):
        super(GenBlastARecord, self).__init__(rangeA, rangeB, float(score))
        self.query = self.rangeA
        self.subject = self.rangeB
        self.cover = float(cover)
        self.aln = int(aln)
        self.rank = int(rank)
        self.hsp_list = hsp_list


def read_genblasta_results(genblasta_output_file):
    try:
        with open(genblasta_output_file, 'r') as f:
            for each_line in f:
                each_line = each_line.strip()
                # record start
                matchObj = re.match("^\/\/.*START.*\/\/$", each_line)
                if matchObj:
                    continue

                # record query name
                matchObj = re.match("\/\/for query:\s+(\S+)\s+\/\/", each_line)
                if matchObj:
                    query_name = matchObj.groups()[0]
                    gba_range_list = []
                    continue

                # subject record
                matchObj = re.match(
                    "^(.*)\|(.*):(\d+)\.\.(\d+)\|(\S+)\|gene cover:(\d+)\((.*)%\)\|score:(.*)\|rank:(\d+)$",
                    each_line)
                if matchObj:
                    if not len(gba_range_list) == 0:
                        gba_range = gba_range_list[-1]
                        gba_range.hsp_list = sorted(
                            gba_range.hsp_list, key=lambda hsp: hsp.query.start)
                        q_range_start = min(
                            [hsp.query.start for hsp in gba_range.hsp_list])
                        q_range_end = max(
                            [hsp.query.end for hsp in gba_range.hsp_list])
                        gba_range.query.range = (q_range_start, q_range_end)

                    q_name, s_name, s_start, s_end, strand, aln, cover, score, rank = matchObj.groups()
                    q_range = ChrLoci(q_name, "+")
                    s_range = ChrLoci(s_name, strand, s_start, s_end)
                    gba_range_list.append(GenBlastARecord(
                        q_range, s_range, score, aln, float(cover) / 100, rank, []))
                    continue

                # hsp record
                matchObj = re.match(
                    "^HSP_ID\[(\d+)\]:\((\d+)-(\d+)\);query:\((\d+)-(\d+)\); pid: (.*)$",
                    each_line)
                if matchObj:
                    hsp_id, s_start, s_end, q_start, q_end, pid = matchObj.groups()
                    q_range = ChrLoci(q_name, "+", q_start, q_end)
                    s_range = ChrLoci(s_name, strand, s_start, s_end)
                    hsp_tmp = BlastHspRecord(q_range, s_range, hsp_id, pid)
                    gba_range_list[-1].hsp_list.append(hsp_tmp)
                    continue

                # record end
                matchObj = re.match("^\/\/.*END.*\/\/$", each_line)
                if matchObj:
                    if not len(gba_range_list) == 0:
                        gba_range = gba_range_list[-1]
                        gba_range.hsp_list = sorted(
                            gba_range.hsp_list, key=lambda hsp: hsp.query.start)
                        q_range_start = min(
                            [hsp.query.start for hsp in gba_range.hsp_list])
                        q_range_end = max(
                            [hsp.query.end for hsp in gba_range.hsp_list])
                        gba_range.query.range = (q_range_start, q_range_end)

                    yield query_name, gba_range_list
                    continue
    except:
        raise ValueError("failed parse %s" % genblasta_output_file)


def gba_record_modify(gba_record, refname_function=None, site_change_function=None):
    if refname_function is not None:
        gba_record.query.chr_id = refname_function(
            gba_record.query.chr_id, "query")
        gba_record.subject.chr_id = refname_function(
            gba_record.subject.chr_id, "subject")
        for hsp in gba_record.hsp_list:
            hsp.query.chr_id = refname_function(hsp.query.chr_id, "query")
            hsp.subject.chr_id = refname_function(
                hsp.subject.chr_id, "subject")

    if site_change_function is not None:
        gba_record.query.range = site_change_function(
            gba_record.query.range, "query")
        gba_record.subject.range = site_change_function(
            gba_record.subject.range, "subject")
        for hsp in gba_record.hsp_list:
            hsp.query.range = site_change_function(hsp.query.range, "query")
            hsp.subject.range = site_change_function(
                hsp.subject.range, "subject")


class GenBlastAJob(object):
    def __init__(self, query_seq_list, subject_seq_list, tmp_work_dir=tmp_work_dir, blast_type='p2n'):
        self.work_dir = tmp_work_dir + "/" + uuid.uuid1().hex
        mkdir(self.work_dir)
        self.data_dir = {
            "query": [self.work_dir + "/query.fa", query_seq_list],
            "subject": [self.work_dir + "/subject.fa", subject_seq_list],
            "output": [self.work_dir + "/output.gba.txt"]
        }

        # rename sequence in query and subject file
        for qs in ["query", "subject"]:
            rename_map = {}
            num = 0
            with open(self.data_dir[qs][0], "w") as f:
                for seq_tmp in self.data_dir[qs][1]:
                    new_name = qs + str(num)
                    rename_map[new_name] = {
                        "old": seq_tmp.seqname, "new": new_name}
                    f.write(">%s\n%s\n" % (new_name, seq_tmp.seq))
                    num = num + 1
            self.data_dir[qs].append(rename_map)

        # set genblasta parameter
        if blast_type == 'p2n':
            self.genblastA_args = {
                "-P": "blast",
                "-pg": "tblastn",
                "-q": self.data_dir['query'][0],
                "-t": self.data_dir['subject'][0],
                "-p": "T",
                "-e": float("1e-2"),
                "-g": "T",
                "-f": "F",
                "-a": 0.5,
                "-d": 100000,
                "-r": 10,
                "-c": 0.5,
                "-s": 0,
                "-o": self.data_dir['output'][0]}
        elif blast_type == 'n2n':
            self.genblastA_args = {
                "-P": "blast",
                "-pg": "blastn",
                "-q": self.data_dir['query'][0],
                "-t": self.data_dir['subject'][0],
                "-p": "F",
                "-e": float("1e-2"),
                "-g": "T",
                "-f": "F",
                "-a": 0.5,
                "-d": 100000,
                "-r": 10,
                "-c": 0.5,
                "-s": 0,
                "-o": self.data_dir['output'][0]}

    def set_para(self, new_para_dict):
        for i in new_para_dict:
            self.genblastA_args[i] = new_para_dict[i]

    def run(self, dry_run=False):
        cmd_string = genblasta_path + " -P %s -pg %s -q %s -t %s -p %s -e %.0e -g %s -f %s -a %.1f -d %d -r %d -c %.1f -s %d -o %s" % (
            self.genblastA_args["-P"],
            self.genblastA_args["-pg"],
            self.genblastA_args["-q"],
            self.genblastA_args["-t"],
            self.genblastA_args["-p"],
            self.genblastA_args["-e"],
            self.genblastA_args["-g"],
            self.genblastA_args["-f"],
            self.genblastA_args["-a"],
            self.genblastA_args["-d"],
            self.genblastA_args["-r"],
            self.genblastA_args["-c"],
            self.genblastA_args["-s"],
            self.genblastA_args["-o"],
        )

        self.gba_cmd = cmd_string

        if not dry_run:
            cmd_run(self.gba_cmd, silence=True)

    def get_old_name(self, new_name, q_or_s):
        if "query" == q_or_s:
            return self.data_dir["query"][2][new_name]['old']
        elif "subject" == q_or_s:
            return self.data_dir["subject"][2][new_name]['old']

    def parse_output(self, rename_block=True):
        for q_name, gba_record_list in read_genblasta_results(self.data_dir['output'][0]):
            if rename_block:
                [gba_record_modify(gba_record, refname_function=self.get_old_name)
                 for gba_record in gba_record_list]
                q_name = self.get_old_name(q_name, "query")
            yield q_name, gba_record_list

    def write_rename_dict(self, output_prefix=None):
        if output_prefix is None:
            output_prefix = self.work_dir + "/seq.rename"

        for qs in ["query", "subject"]:
            output_file = "%s.%s.map" % (output_prefix, qs)
            with open(output_file, 'w') as f:
                for i in self.data_dir[qs]:
                    rename_map_dict = i[2]
                    for j in rename_map_dict:
                        f.write("%s\t%s\n" % (
                            rename_map_dict[j]['new'], rename_map_dict[j]['old']))

    def clean(self):
        shutil.rmtree(self.work_dir)


def genblasta_run(query_file, target_file, output_file, coverage=0.5, rank_num=100, blast_type='wu_blast'):
    if blast_type == 'wu_blast':
        #         genblasta_tmp_dir = "/tmp/genblasta_xyx"
        #         mkdir(genblasta_tmp_dir)
        #
        #         tmp_work_dir = genblasta_tmp_dir + "/" + uuid.uuid1().hex
        #         mkdir(tmp_work_dir)
        #         # tmp_work_dir = os.path.abspath(tmp_work_dir)
        #
        #         blosum62_file = wu_blast_path + "/matrix/aa/BLOSUM62"
        #         tblastn_path = wu_blast_path + "/tblastn"
        #         xdformat_path = wu_blast_path + "/xdformat"
        #         ab_blastall_path = wu_blast_path + "/ab-blastall"
        #
        #         with open(tmp_work_dir + "/cmd.sh", 'w') as f:
        #             f.write("""
        # ln -s %s blosum62
        # ln -s %s tblastn
        # ln -s %s xdformat
        # ln -s %s wu-blastall
        #
        # export PATH=%s:$PATH
        #
        # %s -P wublast -pg tblastn -q %s -t %s -p T -e 1e-2 -g T -f F -a 0.5 -d 100000 -r %d -c %f -s 0 -o %s
        # """ % (blosum62_file, tblastn_path, xdformat_path, ab_blastall_path, tmp_work_dir, genblasta_path, query_file,
        #        target_file, rank_num, coverage, output_file))
        #
        #         cmd_string = "bash cmd.sh"
        #         cmd_run(cmd_string, cwd=tmp_work_dir, silence=True)
        #
        #         rmdir(tmp_work_dir)
        cmd_string = "python %s -o %s -c %f -n %d %s %s" % (
            genblasta_wu_script, output_file, coverage, rank_num, query_file, target_file)
        cmd_run(cmd_string, silence=True)

    else:
        cmd_string = "%s -P blast -pg tblastn -q %s -t %s -p T -e 1e-2 -g T -f F -a 0.5 -d 100000 -r %d -c %f -s 0 -o %s" % (
            genblasta_path, query_file, target_file, rank_num, coverage, output_file)
        cmd_run(cmd_string, silence=True)


def genblasta_to_genewise(query_file, target_file, genblasta_output, work_tmp_dir, flank_range=1000,
                          output_prefix=None, keep=True, genblasta_num=100):
    try:
        cmd_string = "rm -rf core.*"
        cmd_run(cmd_string, retry_max=1, silence=True)
    except:
        pass

    def read_genblasta_results(gen_file):
        with open(gen_file, 'r') as f:
            query_record = {}
            query_record['query'] = ""
            query_record['subject'] = {}
            for each_line in f:
                each_line = each_line.strip()
                # record start
                matchObj = re.match("^\/\/.*START.*\/\/$", each_line)
                if matchObj:
                    continue

                # record query name
                matchObj = re.match("\/\/for query:\s+(\S+)\s+\/\/", each_line)
                if matchObj:
                    query_record['query'] = matchObj.groups()[0]
                    continue

                # subject record
                matchObj = re.match(
                    "^(.*)\|(.*):(\d+)\.\.(\d+)\|(\S+)\|gene cover:(\d+)\((.*)%\)\|score:(.*)\|rank:(\d+)$",
                    each_line)
                if matchObj:
                    q_name, s_name, s_start, s_end, strand, aln, cover, score, rank = matchObj.groups()
                    query_record['subject'][len(query_record['subject'])] = (
                        q_name, s_name, s_start, s_end, strand, aln, cover, score, rank)
                    continue

                # record end
                matchObj = re.match("^\/\/.*END.*\/\/$", each_line)
                if matchObj:
                    yield query_record
                    query_record = {}
                    query_record['query'] = ""
                    query_record['subject'] = {}
                    continue

    try:
        query_seq = Fasta(query_file)
        target_seq = Fasta(target_file)
        gff_dict = {}
        pro_dict = {}
        cds_dict = {}
        failed_job = {}
        failed_id = 0
        for query_record in read_genblasta_results(genblasta_output):
            q_name = query_record['query']
            s_dir = query_record['subject']
            rank_list = []

            for s_id in s_dir:
                if len(rank_list) > genblasta_num:
                    continue

                q_name, s_name, s_start, s_end, strand, aln, cover, score, rank = s_dir[s_id]
                s_start = int(s_start)
                s_end = int(s_end)
                if strand == "+":
                    strand = 1
                else:
                    strand = -1

                uuid_name = uuid.uuid1().hex

                q_file = work_tmp_dir + "/" + uuid_name + ".qfa"
                q_len = len(query_seq[q_name])

                with open(q_file, 'w') as q:
                    q.write(">%s\n%s" %
                            ("query", str(query_seq[q_name]).upper()))

                s_file = work_tmp_dir + "/" + uuid_name + ".sfa"
                s_flank_start = max(min(s_start, s_end) - flank_range, 1)
                s_flank_end = min(max(s_start, s_end) +
                                  flank_range, len(target_seq[s_name]))
                s_seq = str(target_seq[s_name][s_flank_start - 1:s_flank_end])
                s_base = s_flank_start
                if strand == -1:
                    s_seq = reverse_complement(s_seq)
                    s_base = s_flank_end

                with open(s_file, 'w') as s:
                    s.write(">%s\n%s" % ("subject", s_seq.upper()))

                wise_file = work_tmp_dir + "/" + uuid_name + ".gws"

                if rank in rank_list:
                    h_name = str(q_name) + "_" + str(rank) + \
                        "_" + str(rank_list.count(rank) + 1)
                else:
                    h_name = str(q_name) + "_" + str(rank)
                rank_list.append(rank)

                try:
                    genewise_run(q_file, s_file, wise_file, "+")
                    block_dict, gff_output, block_pro_seq, block_cds_seq = genewise_output_parser(
                        wise_file)
                    contig_rec = block_dict_2_gff_rec(
                        block_dict, gff_output, q_len)
                    renamed_contig_rec = gff_rec_obj_rename(contig_rec, hit_name=h_name, q_name=q_name, s_name=s_name,
                                                            s_base=s_base,
                                                            strand=strand)

                    if s_name not in gff_dict:
                        rec = SeqRecord(Seq(""), s_name)
                        gff_dict[s_name] = rec

                    gff_dict[s_name].features.append(
                        renamed_contig_rec.features[0])

                    pro_dict[h_name] = block_pro_seq
                    cds_dict[h_name] = block_cds_seq

                    # clean
                    if not keep:
                        os.remove(s_file)
                        os.remove(q_file)
                        os.remove(wise_file)
                except:
                    failed_id = failed_id + 1
                    failed_job[failed_id] = (
                        q_name, q_file, s_file, wise_file, s_dir[s_id])

        if len(failed_job) > 0:
            return "error"

        gff_file_list = []
        cds_file_list = []
        aa_file_list = []

        if output_prefix is not None:
            cds_file_list.append(output_prefix + ".cds")
            with open(output_prefix + ".cds", 'w') as f:
                for i in cds_dict:
                    f.write(">%s.cds\n%s\n" % (i, cds_dict[i]))

            aa_file_list.append(output_prefix + ".pep")
            with open(output_prefix + ".pep", 'w') as f:
                for i in pro_dict:
                    f.write(">%s.pep\n%s\n" % (i, pro_dict[i]))

            gff_file_list.append(output_prefix + ".gff3")
            rec_obj_list = [gff_dict[i] for i in gff_dict]
            with open(output_prefix + ".gff3", "w") as f:
                GFF.write(rec_obj_list, f)

        return gff_file_list, cds_file_list, aa_file_list

    except:
        return "error"


def gene_to_aa(gene_name):
    return gene_name + ".pep"


def gene_to_cds(gene_name):
    return gene_name + ".cds"


# just work for orthofinder
def query_gene_to_speci(gene_name):
    return gene_name.split("_")[0]


def get_HPR_from_genewise_results(gff_file_list, aa_file_list, cds_file_list, query_file_list, gene_to_aa,
                                  gene_to_cds,
                                  query_gene_to_speci,
                                  subject_species=None):
    """
    get HomoPredictResults class from genblasta_to_genewise output

    seq_parse_way can be xyx or pyfaidx
    """

    query_dict = {}
    for query_file in query_file_list:
        record_dict = read_fasta(query_file)[0]
        for record_id in record_dict:
            record = record_dict[record_id]
            query_dict[record.seqname] = record
    aa_dict = {}
    for aa_file in aa_file_list:
        record_dict = read_fasta(aa_file)[0]
        for record_id in record_dict:
            record = record_dict[record_id]
            aa_dict[record.seqname] = record
    cds_dict = {}
    for cds_file in cds_file_list:
        record_dict = read_fasta(cds_file)[0]
        for record_id in record_dict:
            record = record_dict[record_id]
            cds_dict[record.seqname] = record

    hit_gene_dir = {}
    for gff_file in gff_file_list:
        with open(gff_file, 'r') as in_handle:
            for rec in GFF.parse(in_handle):
                for gene in rec.features:
                    # print(gene.id, gff_file, rec.id)
                    chr_loci_tmp = ChrLoci(rec.id, gene.strand, gene.location.start.position + 1,
                                           gene.location.end.position)
                    hit_gene_tmp = Gene(gene.id, chr_loci=chr_loci_tmp, model_aa_seq=aa_dict[gene_to_aa(gene.id)],
                                        model_cds_seq=cds_dict[gene_to_cds(gene.id)])
                    hit_gene_tmp.gff_obj = gene

                    mRNA_qual = gene.sub_features[0].qualifiers
                    query_gene = mRNA_qual['Target'][0]
                    if not query_gene in hit_gene_dir:
                        hit_gene_dir[query_gene] = []
                    hit_gene_dir[query_gene].append(hit_gene_tmp)

    output_list = []
    for query_id in hit_gene_dir:
        query_gene = Gene(
            query_id, model_aa_seq=query_dict[query_id], species=query_gene_to_speci(query_id))
        output_list.append(
            HomoPredictResults(query_gene, subject_species=subject_species, hit_gene_list=hit_gene_dir[query_id]))

    return output_list


if __name__ == '__main__':
    query_file = "/lustre/home/xuyuxing/Work/tmp/query.faa"
    target_file = "/lustre/home/xuyuxing/Work/tmp/Cuscuta.genome.v1.1.fasta"

    # genblasta object useage

    query_seq_list = (read_fasta_by_faidx(query_file)[
                      i] for i in read_fasta_by_faidx(query_file))
    target_seq_list = (read_fasta_by_faidx(target_file)[
                       i] for i in read_fasta_by_faidx(target_file))

    genblasta_job = GenBlastAJob(
        query_seq_list, target_seq_list, "/lustre/home/xuyuxing/Work/tmp/tmp")
    genblasta_job.run()
    genblasta_out_dir = {i[0]: i[1] for i in genblasta_job.parse_output(True)}

    for query_name in genblasta_out_dir:
        for genblasta_record in genblasta_out_dir[query_name]:
            subject_seq = genblasta_record.subject.get_sequence(target_file)
            print(subject_seq.seqname)
            print(subject_seq.seq)

    # genblasta to genewise
    query_file = "/lustre/home/xuyuxing/Database/Balanophora/Comp_genome/dongming/orthofinder_191218/output/new_run/Cau_gene_loss2/tmp/test/ref_query_19266.fa"
    target_file = "/lustre/home/xuyuxing/Database/Balanophora/Comp_genome/Cau/Cau.genome.fasta"
    genblasta_output = "/lustre/home/xuyuxing/Database/Balanophora/Comp_genome/dongming/orthofinder_191218/output/new_run/Cau_gene_loss2/tmp/test/ref_query_19266.fa.gba"
    work_tmp_dir = "/lustre/home/xuyuxing/Database/Balanophora/Comp_genome/dongming/orthofinder_191218/output/new_run/Cau_gene_loss2/tmp/test/tmp"
    output_prefix = "/lustre/home/xuyuxing/Database/Balanophora/Comp_genome/dongming/orthofinder_191218/output/new_run/Cau_gene_loss2/tmp/test/ref_query_19266"
    flank_range = 1000

    genblasta_run(query_file, target_file, genblasta_output)
    gff_file_list, cds_file_list, aa_file_list = genblasta_to_genewise(query_file, target_file, genblasta_output,
                                                                       work_tmp_dir, flank_range=flank_range,
                                                                       output_prefix=output_prefix, keep=True)

    annotation_output = get_HPR_from_genewise_results(gff_file_list, aa_file_list, cds_file_list, [query_file],
                                                      gene_to_aa, gene_to_cds,
                                                      query_gene_to_speci, subject_species='Cau')

    # debug
    target_file = '/lustre/home/xuyuxing/Work/Gel/Gene_Loss/monocot/gene_loss/Gel.genome.fasta'
    query_file = '/lustre/home/xuyuxing/Work/Gel/Gene_Loss/monocot/gene_loss/Gel_gene_loss/tmp/ref_query_1.fa'
    output_file = '/lustre/home/xuyuxing/Work/Gel/Gene_Loss/monocot/gene_loss/Gel_gene_loss/tmp/ref_query_1.fa.gba'
