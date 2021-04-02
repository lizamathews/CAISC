#! /bin/bash
# This script is used to unarchive .sra files, extract fastq files, 
# split them into 2 separate files, and gzip them. Each line
# represents a different sample.
module load sratoolkit
fastq-dump --split-files --origfmt --gzip SRR1294684.sra
fastq-dump --split-files --origfmt --gzip SRR1294685.sra
fastq-dump --split-files --origfmt --gzip SRR1294686.sra
fastq-dump --split-files --origfmt --gzip SRR1294687.sra
fastq-dump --split-files --origfmt --gzip SRR1294688.sra
fastq-dump --split-files --origfmt --gzip SRR1294689.sra
fastq-dump --split-files --origfmt --gzip SRR1294690.sra
fastq-dump --split-files --origfmt --gzip SRR1294691.sra
fastq-dump --split-files --origfmt --gzip SRR1294692.sra
fastq-dump --split-files --origfmt --gzip SRR1294693.sra
fastq-dump --split-files --origfmt --gzip SRR1294694.sra
fastq-dump --split-files --origfmt --gzip SRR1294695.sra
fastq-dump --split-files --origfmt --gzip SRR1294696.sra
fastq-dump --split-files --origfmt --gzip SRR1294697.sra
fastq-dump --split-files --origfmt --gzip SRR1294698.sra
fastq-dump --split-files --origfmt --gzip SRR1294699.sra