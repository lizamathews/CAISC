# Read in vcf file
rm(list=ls())
library(data.table)
library(tidyr)
library(dplyr)

# Read in your vcf file
#dlist=fread('GSE75688_merged_AllPTS_20mutations.vcf',header=TRUE,skip='#CHROM')
#dlist=fread('GSE75688_merged_AllPTS_30mutations.vcf',header=TRUE,skip='#CHROM')
dlist=fread('GSE75688_merged_AllPTS_40mutations.vcf',header=TRUE,skip='#CHROM')
# Extract only autosome SNPs
#chr_list=sapply(seq(1,22),function(x)paste0('chr',x),simplify = T)
####This line has problem with chr**** set.seed(666);dlist1<-dlist[sample(1:1832506, 200), 1:100]  dlist <- dlist1 %>% filter(`#CHROM`%in%chr_list); 
#dlist2 <- dlist1[dlist1[,`#CHROM`] %in% 1:22, ]  ##dlist <- dlist %>% filter(`#CHROM`%in% 1:2)
#dlist <- dlist1[dlist1[,`#CHROM`] %in% 1:22, ]  ##dlist <- dlist %>% filter(`#CHROM`%in% 1:2)
dlist <- dlist %>% filter(`#CHROM`%in% paste0("chr", 1:22))
Format <- dlist%>%select(FORMAT)
Info <- dlist %>% select(`#CHROM`,POS,REF,ALT)
dlist <- dlist %>% select(-`#CHROM`,-POS,-ID,-REF,-ALT,-QUAL,-FILTER,-INFO,-FORMAT)

# By checking the head of the format, we decided the position of GT (genotype info), AD (allele read depth), and DP (total read depth) separated by ":"
# As shown by the example below, GT position is 1, AD position is 2 and DP position is 3.
# Different version of GATK tool may assign different names of the varaibles. Please refer to the GATK manual for more information.
head(Format)
GT_pos=1
AD_pos=2
DP_pos=3

exinfo_x <- function(info){
  return(sapply(strsplit(info,':'),function(x){ifelse(x[GT_pos]=='./.',NA,as.numeric(strsplit(x[AD_pos],',')[[1]][2]))},simplify=T))
}
exinfo_n <- function(info){
  return(sapply(strsplit(info,':'),function(x){ifelse(x[GT_pos]=='./.',NA,as.numeric(x[DP_pos]))},simplify=T))
}
exinfo_z <- function(info){
  return(sapply(strsplit(info,':'),function(x){ifelse(x[GT_pos]=='./.',NA,ifelse(x[GT_pos]!='0/0',1,0))},simplify=T))
}

# Because the vector is too big. I run it on the cluster instead. All you need is to store the object dlist
X=sapply(dlist,exinfo_x)
N=sapply(dlist,exinfo_n)
Z=sapply(dlist,exinfo_z)
cells<-colnames(X)
cells<-gsub("sampleS", "",unlist(strsplit(cells, "\\."))[1:677*3-2])
colnames(X)<-cells;colnames(N)<-cells;colnames(Z)<-cells;
# Quick fitering
# If 99% of the cells did not detect any read, we remove such entry.
thres=floor(0.01*ncol(X))
sel=rowSums(!is.na(Z))>thres
X=X[sel,]
N=N[sel,]
Z=Z[sel,]
#save(Info,X,N,Z,file='DENDRO_input_20mutations.rda')
#save(Info,X,N,Z,file='DENDRO_input_30mutations.rda')
save(Info,X,N,Z,file='DENDRO_input_40mutations.rda')
