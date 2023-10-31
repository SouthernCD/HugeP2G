# HugeP2G
An integration tool for aligning large numbers of protein sequences to the genome (use genBlastA and GeneWise2)

GeneWise2 gives very fine amino acid to genome alignment results through dynamic programming algorithms, but it is very inefficient. genBlastA allows you to align the amino acid sequences to the approximate region of the genome first, and then go through GeneWise2 to give a precise alignment. HugeP2G provides an automated process for this.

## Installation
You need genBlastA and GeneWise2 to run HugeP2G, in which genBlastA is included in the HugeP2G package, and GeneWise2 can be downloaded from [here](https://github.com/SouthernCD/genewise2_docker).

Install HugeP2G from PyPI:
```
pip install HugeP2G
```

## Usage

```
usage: HugeP2G [-h] [-s SKIP_RANGE_FILE] [-d WORK_DIR] [-t NUM_THREADS] [-c GENE_COVERAGE] [-n GENBLASTA_HIT_NUM] [-sc SKIP_COVERAGE] [-split SEQ_NUM_IN_SUBDIR] [-r] query_protein_table target_genome_fasta

Aligning a large number of protein sequences to a genome (use genblasta and genewise)

positional arguments:
  query_protein_table   Path of query genome table in tsv format, must have column "sp_id" and "pt_file", "sp_id" is the species id, "pt_file" is the path of protein fasta file
  target_genome_fasta   Path of target genome fasta file

optional arguments:
  -h, --help            show this help message and exit
  -s SKIP_RANGE_FILE, --skip_range_file SKIP_RANGE_FILE
                        Path of skip_range_file, tsv file, should have column name "chr", "start", and "end", the range of the genome that need to be skipped
  -d WORK_DIR, --work_dir WORK_DIR
                        Path of work dir (default as ./hugep2g_out)
  -t NUM_THREADS, --num_threads NUM_THREADS
                        threads number (default as 56)
  -c GENE_COVERAGE, --gene_coverage GENE_COVERAGE
                        gene coverage (default as 0.2)
  -n GENBLASTA_HIT_NUM, --genblasta_hit_num GENBLASTA_HIT_NUM
                        genblasta hit num (default as 50)
  -sc SKIP_COVERAGE, --skip_coverage SKIP_COVERAGE
                        annotated_coverage (default as 0.8)
  -split SEQ_NUM_IN_SUBDIR, --seq_num_in_subdir SEQ_NUM_IN_SUBDIR
                        split big fasta file to run (default as 1000)
  -r, --force_redo      force redo all job
```

## Example
```
HugeP2G -t 80 query.tsv genome.fasta
```