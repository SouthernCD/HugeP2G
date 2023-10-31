import argparse
from hugep2g.pipelines import build_hugep2g


def main():

    # argument parse
    parser = argparse.ArgumentParser(
        prog='HugeP2G', description='Aligning a large number of protein sequences to a genome (use genblasta and genewise)'
    )

    parser.add_argument("query_protein_table",
                        help="Path of query genome table in tsv format, must have column \"sp_id\" and \"pt_file\", \"sp_id\" is the species id, \"pt_file\" is the path of protein fasta file",
                        type=str)
    parser.add_argument("target_genome_fasta",
                        help="Path of target genome fasta file", type=str)
    parser.add_argument("-d", "--work_dir", help="Path of work dir (default as ./hugep2g_out)",
                        default="./hugep2g_out", type=str)
    parser.add_argument(
        "-t", "--num_threads", help="threads number (default as 56)", default=56, type=int)
    parser.add_argument("-c", "--gene_coverage",
                        help="min gene coverage (default as 0.2)", default=0.2, type=float)
    parser.add_argument("-n", "--genblasta_hit_num",
                        help="max genblasta hit num (default as 50)", default=50, type=int)
    parser.add_argument("-s", "--skip_range_file",
                        help="Path of skip_range_file, tsv file, should have column name \"chr\", \"start\", and \"end\", the range of the genome that need to be skipped", type=str)
    parser.add_argument("-sc", "--skip_coverage", help="Hits that exceed this ratio with respect to the area coverage given in skip_range_file will be ignored (default as 0.8)", default=0.8,
                        type=float)
    parser.add_argument("-split", "--seq_num_in_subdir", help="split big fasta file to run (default as 5000)", default=5000,
                        type=int)
    parser.add_argument("-r", "--force_redo",
                        help="force redo all job", action='store_true')

    args = parser.parse_args()

    build_hugep2g(args)

if __name__ == '__main__':
    main()