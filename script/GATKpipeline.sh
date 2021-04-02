#! /bin/bash
#SBATCH --job-name="subread"
#SBATCH --mail-type=BEGIN,END
#SBATCH --mem=50g
#SBATCH --cpus-per-task=50
#SBATCH --time=10:00:00
#$ -j y
#$ -cwd
#$ -pe orte 35
#$ -q main.q@miner4

### Modify file paths as necessary

# path where packages are
data=$3

# Run this line only once, and then uncomment when you're done
data/GATKpackages/STAR-STAR_2.4.0k/bin/Linux_x86_64_static/STAR --runMode genomeGenerate --genomeDir data/STARGRCh37index --genomeFastaFiles data/ucsc_hg19/ucsc.hg19.fasta --runThreadN 10
module load java/1.7.0_25 # run this line if on cluster

# Run from command line: bash GATKpipeline.sh sample fastq_path data

# path where fastq files are
fastq_path=$2

##Alignment jobs were executed as follows:
sample=$1
echo ${sample}

mkdir ${sample}

gunzip fastq_path/${sample}_R1.fastq.gz
gunzip fastq_path/${sample}_R2.fastq.gz
cd ${sample}

runDir=${sample}_1pass
mkdir $runDir
cd $runDir
data/GATKpackages/STAR-STAR_2.4.0k/bin/Linux_x86_64_static/STAR --genomeDir data/STARGRCh37index --readFilesIn fastq_path/${sample}_R1.fastq fastq_path/${sample}_R2.fastq --runThreadN 20
cd ..

#For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass:
genomeDirS=${sample}_2pass
mkdir $genomeDirS
#cd $genomeDirS
#cd ..

data/GATKpackages/STAR-STAR_2.4.0k/bin/Linux_x86_64_static/STAR --runMode genomeGenerate --genomeDir $genomeDirS --genomeFastaFiles data/ucsc_hg19/ucsc.hg19.fasta --sjdbFileChrStartEnd ${runDir}/SJ.out.tab --sjdbOverhang 75 --runThreadN 20

#The resulting index is then used to produce the final alignments as follows:
data/GATKpackages/STAR-STAR_2.4.0k/bin/Linux_x86_64_static/STAR --genomeDir $genomeDirS --readFilesIn  fastq_path/${sample}_R1.fastq fastq_path/${sample}_R2.fastq --runThreadN 20

###The above step produces a SAM file, which we then put through the usual Picard processing steps: adding read group information, sorting, marking duplicates and indexing.
java -jar data/GATKpackages/GATK3.3/picard-tools-1.119/AddOrReplaceReadGroups.jar I=Aligned.out.sam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sampleS${sample}
java -jar data/GATKpackages/GATK3.3/picard-tools-1.119/MarkDuplicates.jar I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

###Splicing N trim and quality conversion
java -jar data/GATKpackages/GATK3.3/GenomeAnalysisTK.jar -T SplitNCigarReads -R data/ucsc_hg19/ucsc.hg19.fasta -I dedupped.bam -o ${sample}_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

rm Aligned.out.sam
rm rg_added_sorted.bam
rm rg_added_sorted.bam.bai
rm dedupped.bam
rm dedupped.bai
rm $genomeDirS/SA
rm $genomeDirS/Genome
rm $genomeDirS/SAindex
rm $runDir/Aligned.out.sam

#Variant calling
java -jar data/GATKpackages/GATK3.3/GenomeAnalysisTK.jar -T HaplotypeCaller -R data/GATKpackages/ucsc_hg19/ucsc.hg19.fasta -I ${sample}_split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o ${sample}_output.vcf
#Variant filtering
java -jar data/GATKpackages/GATK3.3/GenomeAnalysisTK.jar -T VariantFiltration -R data/GATKpackages/ucsc_hg19/ucsc.hg19.fasta -V ${sample}_output.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${sample}_output_final.vcf
cd ..
gzip fastq_path/${sample}_R1.fastq
gzip fastq_path/${sample}_R2.fastq




