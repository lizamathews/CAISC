#!/bin/bash
#SBATCH --job-name="merge_vcfs"
#SBATCH --mem=25g
#SBATCH --time=2:00:00
#SBATCH --partition="quick"

# Run from command line if not using cluster: bash merge_vcf.sh file_path output_vcf
# Run from command line if using cluster: sbatch merge_vcf.sh file_path output_vcf

# Merge VCF Files using List of VCF Files

cd $1/results

#module load java/1.7.0_25 # run this line if on cluster

java -jar $1/GATK3.3/GenomeAnalysisTK.jar -T CombineVariants -R $1/ucsc_hg19/ucsc.hg19.fasta --variant $1/results/filelist.list -o $2 -genotypeMergeOptions UNIQUIFY
