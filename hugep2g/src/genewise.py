from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from hugep2g.config import genewise_path
from toolbiox.lib.common.os import cmd_run
import re


def genewise_run(query_file, subject_file, genewise_out, strand, pseudo=False, q_start=None, q_end=None):
    if pseudo:
        if strand == "+":
            cmd_string = genewise_path + " -pretty -tfor -genes -gff -pseudo -pep -nosplice_gtag -init global -s %d -t %d %s %s > %s" % (
                q_start, q_end, query_file, subject_file, genewise_out)
        else:
            cmd_string = genewise_path + " -pretty -trev -genes -gff -pseudo -pep -nosplice_gtag -init global -s %d -t %d %s %s > %s" % (
                q_start, q_end, query_file, subject_file, genewise_out)
    else:
        if strand == "+":
            cmd_string = genewise_path + " -pretty -tfor -genes -gff -pseudo -pep %s %s > %s" % (
                query_file, subject_file, genewise_out)
        else:
            cmd_string = genewise_path + " -pretty -trev -genes -gff -pseudo -pep %s %s > %s" % (
                query_file, subject_file, genewise_out)

    cmd_run(cmd_string, silence=True)

    return query_file, subject_file, genewise_out


def parse_align_sub(align_sub):
    p0 = re.compile('^$')
    p00 = re.compile('^\s+$')
    p_1 = re.compile('(...................\d) (.*)$')
    p_23 = re.compile('(                    ) (.*)$')
    p_wise_mark = re.compile('^.*.[^\d](\d+)')
    True_block = []
    num = 0
    flag1 = 0
    flag2 = 0
    flag3 = 0
    flag4 = 0
    flag5 = 0
    flag6 = 0
    line1 = []
    line2 = []
    line3 = []
    line4 = []
    line5 = []
    line6 = []
    for block in align_sub:
        block = re.sub(p0, '', block)
        block = re.sub(p00, '', block)
        lines = block.split('\n')
        block_start = 0
        for line in lines:
            if re.match('Score (.*) bits over entire alignment', line):
                genewise_score = float(re.findall(
                    'Score (.*) bits over entire alignment', line)[0])
            elif flag1 == 0:
                match = p_1.match(line)
                if match:
                    query_name = re.findall('^(\S+)\s+.*', line)[0]
                    flag1 = 1
                    wise_mark = match.group(1)
                    wise_str = match.group(2)
                    match2 = p_wise_mark.match(wise_mark)
                    query_block_start = match2.group(1)
                    line1 = wise_str
                    continue
            elif flag1 == 1 and flag2 == 0:
                match = p_23.match(line)
                if match:
                    flag2 = 1
                    wise_str = match.group(2)
                    line2 = wise_str
                    continue
            elif flag1 == 1 and flag2 == 1 and flag3 == 0:
                match = p_23.match(line)
                if match:
                    flag3 = 1
                    wise_str = match.group(2)
                    line3 = wise_str
                    continue
            elif flag1 == 1 and flag2 == 1 and flag3 == 1 and flag4 == 0:
                match = p_1.match(line)
                if match:
                    subject_name = re.findall('^(\S+)\s+.*', line)[0]
                    flag4 = 1
                    wise_mark = match.group(1)
                    wise_str = match.group(2)
                    match2 = p_wise_mark.match(wise_mark)
                    block_start = match2.group(1)
                    line4 = wise_str
                    continue
            elif flag1 == 1 and flag2 == 1 and flag3 == 1 and flag4 == 1 and flag5 == 0:
                match = p_23.match(line)
                if match:
                    flag5 = 1
                    wise_str = match.group(2)
                    line5 = wise_str
                    continue
            elif flag1 == 1 and flag2 == 1 and flag3 == 1 and flag4 == 1 and flag5 == 1 and flag6 == 0:
                match = p_23.match(line)
                if match:
                    flag6 = 1
                    wise_str = match.group(2)
                    line6 = wise_str
                    num = num + 1
                    True_block.append([num, block_start, query_block_start, [
                                      line1, line2, line3, line4, line5, line6]])
                    flag1 = 0
                    flag2 = 0
                    flag3 = 0
                    flag4 = 0
                    flag5 = 0
                    flag6 = 0
                    line1 = []
                    line2 = []
                    line3 = []
                    line4 = []
                    line5 = []
                    line6 = []
                    continue
    return True_block, query_name, subject_name, genewise_score


def True_block_merge(True_block):
    """
    :return: subject_start, query_start, merged_block
    """
    merged_block = ["", "", "", "", "", ""]
    for i in True_block:
        for j in range(0, 6):
            merged_block[j] = merged_block[j] + i[3][j]

    merged_block_clean = ["", "", "", "", "", ""]
    for i in range(0, len(merged_block[0])):
        if set([merged_block[j][i] for j in range(0, 6)]) == set(" "):
            continue
        else:
            for j in range(0, 6):
                merged_block_clean[j] = merged_block_clean[j] + \
                    merged_block[j][i]

    return True_block[0][1], True_block[0][2], merged_block_clean


def cds_parser(cds_block, position=1):
    line1, line2, line3, line4, line5, line6 = cds_block
    num = -1
    stop = []
    frame = []
    cds = ""
    protein = ""
    begin = position
    position = position - 1
    for i in line3:
        num = num + 1
        stop_flag = 0
        pseudo_flag = 0
        # print "\n"
        # print position
        # print i
        # print line4[num]+line5[num]+line6[num]
        if i == "-":  # GAP
            continue
        elif i == " " and line4[num] == " ":  # Null
            continue
        elif i == " " and (not line4[num] == " "):  # frame_intron
            position = position + 1
            cds = cds + line4[num]
        elif i == "X":  # stop code
            stop.append(position)
            position = position + 3
            cds = cds + line4[num] + line5[num] + line6[num]
            protein = protein + "X"
        elif i == "!":  # frame_pseudo
            frame.append(position)
            protein = protein + "X"
            cds = cds + "NNN"
            position = position + int(line4[num])
        else:  # normal
            position = position + 3
            protein = protein + i
            cds = cds + line4[num] + line5[num] + line6[num]

    return cds, protein, frame, stop, position - begin + 1


def intron_parser(merged_block_clean, block_subject_start, query_start):
    p_intron = re.compile('<(\d)-*\[(\d+) *\: *(\d+)\]-*(\d)>')
    line1, line2, line3, line4, line5, line6 = merged_block_clean

    if len(list(re.finditer(p_intron, line5))) == 0:
        block_dict = OrderedDict()
        last_intron_end = -1
        last_intron_end_site = block_subject_start - 1
        last_intron_phase = 0
        num = 0

        # make cds
        last_cds = ["", "", "", "", "", ""]
        num = num + 1
        for i in range(last_intron_end + 1, len(line5)):
            for j in range(0, 6):
                last_cds[j] = last_cds[j] + merged_block_clean[j][i]
        block_dict[num] = {}
        block_dict[num]["block"] = last_cds
        block_dict[num]['type'] = 'cds'
        block_dict[num]['phase'] = last_intron_phase
        block_dict[num]['start'] = last_intron_end_site + 1
        block_dict[num]['cds'], block_dict[num]['pro'], block_dict[num]['frame'], block_dict[num]['stop'], \
            block_dict[num]['true_length'] = cds_parser(
            block_dict[num]["block"], block_dict[num]['start'])

        whole_length = sum([block_dict[i]['true_length'] for i in block_dict])

        block_subject_end = block_subject_start + whole_length - 1
        block_dict[num]['end'] = block_subject_end

        # add query range
        last_q_end = query_start - 1
        for i in block_dict:
            if block_dict[i]['type'] == 'cds':
                block_dict[i]['q_start'] = last_q_end + 1
                block_dict[i]['q_end'] = block_dict[i]['q_start'] + \
                    len(block_dict[i]['block'][0].replace('-', '')) - 1
                last_q_end = block_dict[i]['q_end']

        return block_dict
    else:
        block_dict = OrderedDict()
        last_intron_end = -1
        last_intron_end_site = block_subject_start - 1
        last_intron_phase = 0
        num = 0
        for intron_match in re.finditer(p_intron, line5):
            phase, i_start, i_end, phase = intron_match.groups()
            i_start = int(i_start)
            i_end = int(i_end)
            phase = int(phase)

            # make cds
            last_cds = ["", "", "", "", "", ""]
            num = num + 1
            for i in range(last_intron_end + 1, intron_match.start()):
                for j in range(0, 6):
                    last_cds[j] = last_cds[j] + merged_block_clean[j][i]
            block_dict[num] = {}
            block_dict[num]["block"] = last_cds
            block_dict[num]['type'] = 'cds'
            block_dict[num]['phase'] = last_intron_phase
            block_dict[num]['start'] = last_intron_end_site + 1
            block_dict[num]['end'] = int(i_start) - 1
            block_dict[num]['cds'], block_dict[num]['pro'], block_dict[num]['frame'], block_dict[num]['stop'], \
                block_dict[num]['true_length'] = cds_parser(
                block_dict[num]["block"], block_dict[num]['start'])

            last_intron_end = intron_match.end() - 1
            last_intron_end_site = i_end

            # make intron
            intron_block = ["", "", "", "", "", ""]
            for i in range(intron_match.start(), intron_match.end() + 1):
                for j in range(0, 6):
                    intron_block[j] = intron_block[j] + \
                        merged_block_clean[j][i]

            num = num + 1
            block_dict[num] = {}
            block_dict[num]["block"] = intron_block
            block_dict[num]['type'] = 'intron'
            block_dict[num]['phase'] = 0
            last_intron_phase = (3 - phase) % 3
            # last_intron_phase = phase
            block_dict[num]['start'] = i_start
            block_dict[num]['end'] = i_end
            block_dict[num]['cds'] = ''
            block_dict[num]['true_length'] = block_dict[num]['end'] - \
                block_dict[num]['start'] + 1
            match_obj = re.findall(
                r'(\S):(\S)\[(...)\]', block_dict[num]['block'][2])
            if len(match_obj) == 0 and block_dict[num]['phase'] != 0:
                raise ValueError('phase not as 0, but failed to found aa')
            if len(match_obj) != 0:
                block_dict[num]['pro'] = match_obj[0][1]
            else:
                block_dict[num]['pro'] = ""

        # make tail cds
        last_cds = ["", "", "", "", "", ""]
        num = num + 1
        for i in range(last_intron_end + 1, len(line5)):
            for j in range(0, 6):
                last_cds[j] = last_cds[j] + merged_block_clean[j][i]
        block_dict[num] = {}
        block_dict[num]["block"] = last_cds
        block_dict[num]['type'] = 'cds'
        block_dict[num]['phase'] = last_intron_phase
        block_dict[num]['start'] = last_intron_end_site + 1

        block_dict[num]['cds'], block_dict[num]['pro'], block_dict[num]['frame'], block_dict[num]['stop'], \
            block_dict[num][
            'true_length'] = cds_parser(block_dict[num]["block"], block_dict[num]['start'])

        whole_length = sum([block_dict[i]['true_length'] for i in block_dict])

        block_subject_end = block_subject_start + whole_length - 1
        block_dict[num]['end'] = block_subject_end

        # add query range
        last_q_end = query_start - 1
        for i in block_dict:
            if block_dict[i]['type'] == 'cds':
                block_dict[i]['q_start'] = last_q_end + 1
                block_dict[i]['q_end'] = block_dict[i]['q_start'] + \
                    len(block_dict[i]['block'][0].replace('-', '')) - 1
                last_q_end = block_dict[i]['q_end']

        return block_dict


def genewise_output_parser(genewise_out):
    with open(genewise_out, 'r') as f:
        all_text = f.read()
        info = all_text.split('//')
        # print(len(info))
        # if not len(info) == 5:
        #     raise ValueError("Error: Bad genewise output file\n")

    # align_sub
    align_sub = info[0].split('\n\n')
    True_block, query_name, subject_name, genewise_score = parse_align_sub(
        align_sub)

    # get query start end
    query_start = int(True_block[0][2])
    query_end_temp = int(True_block[-1][2])
    last_block_query = True_block[-1][3][0]
    last_block_query = re.sub(r'\s+', '', last_block_query)
    last_block_query = re.sub(r'-', '', last_block_query)
    query_end = query_end_temp + len(last_block_query) - 1

    # merge and parse block
    subject_start, query_start_block, merged_block = True_block_merge(
        True_block)

    if int(query_start) != int(query_start_block):
        raise ValueError("query start bug")

    block_dict = intron_parser(merged_block, int(subject_start), query_start)
    block_pro_seq = "".join([block_dict[i]['pro'] for i in block_dict]).upper()
    block_cds_seq = "".join([block_dict[i]['cds'] for i in block_dict]).upper()

    # protein_sub
    if len(info[2].split('>')) > 1:
        protein_seq = "".join(info[2].split('>')[1].split('\n')[1:])

    # gff_sub
    gff_sub_cds_dict = {}
    num = 0
    gff_output = []
    if len(info) > 4:
        for i in info[3].split('\n'):
            if i == '':
                continue
            gff_contig, gff_predict, gff_type, gff_start, gff_end, gff_score, gff_strand, gff_phase, gff_attr = i.split(
                '\t')
            if gff_type == 'cds' or gff_type == 'intron':
                num = num + 1
                gff_sub_cds_dict[num] = {}
                gff_sub_cds_dict[num]['type'] = gff_type
                if gff_phase == '.':
                    gff_phase = 0
                gff_sub_cds_dict[num]['phase'] = int(gff_phase)
                gff_sub_cds_dict[num]['start'] = int(gff_start)
                gff_sub_cds_dict[num]['end'] = int(gff_end)
            elif gff_type == 'match':
                gff_output = [gff_contig, float(gff_score), gff_attr, int(gff_start), int(gff_end), query_start,
                              query_end]
    else:
        info_list = info[1].split('\n')
        for i in range(0, len(info_list)):
            string_tmp = info_list[i]
            if string_tmp == '':
                continue
            gff_info = string_tmp.split()
            if gff_info[0] == 'Gene' and len(gff_info) < 3:
                continue
            elif gff_info[0] == 'Gene' and len(gff_info) >= 3:
                gff_output = [subject_name, genewise_score, None, int(gff_info[1]), int(gff_info[2]), query_start,
                              query_end]
            elif len(gff_info) == 5:
                num = num + 1
                gff_sub_cds_dict[num] = {}
                gff_sub_cds_dict[num]['type'] = 'cds'
                gff_sub_cds_dict[num]['phase'] = (3 - int(gff_info[4])) % 3
                gff_sub_cds_dict[num]['start'] = int(gff_info[1])
                gff_sub_cds_dict[num]['end'] = int(gff_info[2])

                if info_list[i + 1] == '':
                    break

                num = num + 1
                gff_sub_cds_dict[num] = {}
                gff_sub_cds_dict[num]['type'] = 'intron'
                gff_sub_cds_dict[num]['phase'] = 0
                gff_sub_cds_dict[num]['start'] = int(gff_info[2]) + 1
                gff_sub_cds_dict[num]['end'] = int(
                    info_list[i + 1].split()[1]) - 1

    # block_check
    # phase have bug in genewise gff sub output, so we can't check it
    # for rank_id in block_dict:
    #     for check_type in ['start', 'end', 'type']:
    #         if not block_dict[rank_id][check_type] == gff_sub_cds_dict[rank_id][check_type]:
    #             raise ValueError("block check failed %d: %s" % (rank_id, check_type))
    # if not block_pro_seq == protein_seq:
    #     raise ValueError("block check failed protein seq %s" % genewise_out)

    return block_dict, gff_output, block_pro_seq, block_cds_seq


def block_dict_2_gff_rec(block_dict, gff_output, query_length):
    subject_id, score, hit_id, start, end, query_start, query_end = gff_output
    # make contig rec
    contig_rec = SeqRecord(Seq(""), "subject")

    # make gene rec
    qualifiers = {
        "source": "genewise",
        "score": score,
        "ID": "genewise_hit",
    }
    top_feature = SeqFeature(FeatureLocation(
        start - 1, end), type="gene", strand=1, qualifiers=qualifiers)

    # make mRNA rec
    mRNA_feature = SeqFeature(FeatureLocation(start - 1, end), strand=1)
    qualifiers = {
        "source": "genewise",
        "score": score,
        "ID": "genewise_hit.mRNA",
        "Target": "query",
        "Target_Start": query_start,
        "Target_End": query_end,
        "Target_Length": query_length
    }

    mRNA_feature.qualifiers = qualifiers
    mRNA_feature.sub_features = []
    mRNA_feature.type = "mRNA"
    top_feature.sub_features = [mRNA_feature]

    # make cds rec
    cds_num = 1
    pseudo_flag_gene = 0
    for rank_id in block_dict:
        if block_dict[rank_id]['type'] == 'cds':
            start = block_dict[rank_id]['start']
            end = block_dict[rank_id]['end']

            feature_now = SeqFeature(FeatureLocation(start - 1, end), strand=1)

            feature_now.type = 'CDS'
            qualifiers = {
                "source": "genewise",
                "score": 0.0,
                "phase": block_dict[rank_id]['phase'],
                "ID": "genewise_hit.cds" + str(cds_num),
                "Target": "query",
                "Target_Start": block_dict[rank_id]['q_start'],
                "Target_End": block_dict[rank_id]['q_end'],
                "Target_Length": query_length
            }
            feature_now.qualifiers = qualifiers

            if len(block_dict[rank_id]['frame']):
                feature_now.qualifiers["Note_frame_shift"] = []
                for j in block_dict[rank_id]['frame']:
                    feature_now.qualifiers["Note_frame_shift"].append(j)
                pseudo_flag_gene = 1
            if len(block_dict[rank_id]['stop']):
                feature_now.qualifiers["Note_stop_code"] = []
                for k in block_dict[rank_id]['stop']:
                    feature_now.qualifiers["Note_stop_code"].append(k)
                pseudo_flag_gene = 1

            top_feature.sub_features[0].sub_features.append(feature_now)

            cds_num = cds_num + 1

    # add pseudo flag
    if pseudo_flag_gene:
        top_feature.sub_features[0].qualifiers["Note"] = "pseudogene"

    contig_rec.features = [top_feature]

    return contig_rec


def gff_rec_obj_rename(gff_rec_object, hit_name="genewise_hit", q_name="query", s_name="subject", s_base=1, strand=1):
    gff_rec_object.id = s_name
    for gene in gff_rec_object.features:
        if strand == 1:
            abs_start = s_base + gene.location.start + 1 - 1
            abs_end = s_base + gene.location.end - 1
        else:
            abs_end = s_base - (gene.location.start + 1) + 1
            abs_start = s_base - gene.location.end + 1
        gene.location = FeatureLocation(abs_start - 1, abs_end)
        gene.strand = strand
        gene.id = hit_name
        gene.qualifiers['ID'] = hit_name
        for mRNA in gene.sub_features:
            if strand == 1:
                abs_start = s_base + mRNA.location.start + 1 - 1
                abs_end = s_base + mRNA.location.end - 1
            else:
                abs_end = s_base - (mRNA.location.start + 1) + 1
                abs_start = s_base - mRNA.location.end + 1
            mRNA.location = FeatureLocation(abs_start - 1, abs_end)
            mRNA.strand = strand
            mRNA.id = hit_name + ".mRNA"
            mRNA.qualifiers['Target'] = q_name
            mRNA.qualifiers['ID'] = mRNA.id
            cds_num = 0
            for cds in mRNA.sub_features:
                cds_num = cds_num + 1
                if strand == 1:
                    abs_start = s_base + cds.location.start + 1 - 1
                    abs_end = s_base + cds.location.end - 1
                else:
                    abs_end = s_base - (cds.location.start + 1) + 1
                    abs_start = s_base - cds.location.end + 1
                cds.location = FeatureLocation(abs_start - 1, abs_end)
                cds.strand = strand
                cds.id = hit_name + ".cds" + str(cds_num)
                cds.qualifiers['ID'] = cds.id
                cds.qualifiers['Target'] = q_name
                for pseudo_type in ["Note_stop_code", "Note_frame_shift"]:
                    if pseudo_type in cds.qualifiers:
                        tmp_new = []
                        for i in cds.qualifiers[pseudo_type]:
                            if strand == 1:
                                abs_site = s_base + int(i) - 1
                            else:
                                abs_site = s_base - int(i) + 1
                            tmp_new.append(abs_site)
                        cds.qualifiers[pseudo_type] = tmp_new

    return gff_rec_object


if __name__ == '__main__':

    subject_file = '/lustre/home/xuyuxing/Work/tmp/C000N_14332_25011.fa'
    query_file = '/lustre/home/xuyuxing/Work/tmp/C000N0003E0.1.aa'
    genewise_out = '/lustre/home/xuyuxing/Work/tmp/C000N0003E0.gws'
    q_name = "C000N0003E0"
    s_name = "C000N"
    h_name = "genewise_hit"
    s_base = 14332
    strand = "+"
    output_prefix = '/lustre/home/xuyuxing/Work/tmp/C000N0003E0.genewise'

    """
    subject_file = '/lustre/home/xuyuxing/Work/tmp/tmp/145ba845d80f4f2fafe174b70675df32.fr.sfa'
    query_file = '/lustre/home/xuyuxing/Work/tmp/tmp/145ba845d80f4f2fafe174b70675df32.qfa'
    genewise_out = '/lustre/home/xuyuxing/Work/tmp/tmp/145ba845d80f4f2fafe174b70675df32.fr.gws'
    output_prefix = '/lustre/home/xuyuxing/Work/tmp/tmp/145ba845d80f4f2fafe174b70675df32.fr.gwse'
    s_base = 1
    """

    genewise_run(query_file, subject_file, genewise_out, strand)
    block_dict, gff_output, block_pro_seq, block_cds_seq = genewise_output_parser(
        genewise_out)

    record_dict = SeqIO.to_dict(SeqIO.parse(query_file, "fasta"))
    query_length = len(str(record_dict[list(record_dict.keys())[0]].seq))
    contig_rec = block_dict_2_gff_rec(block_dict, gff_output, query_length)

    renamed_contig_rec = gff_rec_obj_rename(contig_rec, hit_name=h_name, q_name=q_name, s_name=s_name, s_base=s_base,
                                            strand=1)

    if not output_prefix is None:
        with open(output_prefix + ".cds", 'w') as f:
            f.write(">%s.cds\n%s" % (h_name, block_cds_seq))

        with open(output_prefix + ".pep", 'w') as f:
            f.write(">%s.pep\n%s" % (h_name, block_pro_seq))

        with open(output_prefix + ".gff3", "w") as f:
            GFF.write([renamed_contig_rec], f)
