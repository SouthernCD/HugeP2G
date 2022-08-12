import os
from toolbiox.lib.common.os import cmd_run
from hugep2g.config import project_dir
import argparse


def args_parser():

    # argument parse
    parser = argparse.ArgumentParser(
        prog='genblasta_wu', description='Running GenBlastA with Wu-Blast)'
    )

    parser.add_argument("query_fasta",
                        help="Path of query genome table in excel format (.xlsx), must have column \"id\" and \"pt_file\"",
                        type=str)
    parser.add_argument("target_fasta",
                        help="Path of target genome fasta file", type=str)
    parser.add_argument("-o", "--output_file", help="Path of output file (default as genblasta.txt)",
                        default="genblasta.txt", type=str)
    parser.add_argument("-c", "--gene_coverage",
                        help="gene coverage (default as 0.2)", default=0.2, type=float)
    parser.add_argument("-n", "--genblasta_hit_num",
                        help="genblasta hit num (default as 50)", default=50, type=int)

    args = parser.parse_args()

    args.query_fasta = os.path.realpath(args.query_fasta)
    args.output_file = os.path.realpath(args.output_file)
    args.target_fasta = os.path.realpath(args.target_fasta)

    return args


def main():

    args = args_parser()

    genblasta_dir = os.path.join(project_dir, 'dep', 'genblasta')
    os.environ['PATH'] += os.pathsep + genblasta_dir

    cmd_string = "genblasta_v1.0.4_linux_x86_64 -P wublast -pg tblastn -q %s -t %s -p T -e 1e-2 -g T -f F -a 0.5 -d 100000 -r %d -c %f -s 0 -o %s" % (
        args.query_fasta, args.target_fasta, args.genblasta_hit_num, args.gene_coverage, args.output_file)

    cmd_run(cmd_string, cwd=genblasta_dir, retry_max=1, silence=True)


if __name__ == '__main__':
    main()

    """
    #! /bin/bash

    query=$1
    target=$2
    output=$3
    rank_num=$4
    coverage=$5

    query=`realpath $query`
    target=`realpath $target`
    output=`realpath $output`


    export PATH=/lustre/home/xuyuxing/Program/genblast/genblasta_wu:$PATH

    cd /lustre/home/xuyuxing/Program/genblast/genblasta_wu

    /lustre/home/xuyuxing/Program/genblast/genblasta_v1.0.4_linux_x86_64 -P wublast -pg tblastn -q $query -t $target -p T -e 1e-2 -g T -f F -a 0.5 -d 100000 -r $rank_num -c $coverage -s 0 -o $output
    """
