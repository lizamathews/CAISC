#' Create SNV Distance Matrix input from VCF File
#'
#' This function finds chromosome information, alternative allele read counts (X),
#'  total read counts (N), and mutation profile matrix (Z) from VCF file of mutations
#'
#' @param filevcf Path to input VCF file
#' @return List with Mutation Information (Info, X, N, Z)
vcf_to_SNVinput <- function(filevcf){

  # Read in your vcf file
  dlist=data.table::fread(filevcf,header=TRUE,skip='#CHROM')
  # Extract only autosome SNPs
  dlist <- dlist %>% dplyr::filter(`#CHROM`%in% paste0("chr", 1:22))
  Format <- dlist %>% dplyr::select(FORMAT)
  Info <- dlist %>% dplyr::select(`#CHROM`,POS,REF,ALT)
  dlist <- dlist %>% dplyr::select(-`#CHROM`,-POS,-ID,-REF,-ALT,-QUAL,-FILTER,-INFO,-FORMAT)

  # By checking the head of the format, we decided the position of GT (genotype info), AD (allele read depth), and DP (total read depth) separated by ":"
  # As shown by the example below, GT position is 1, AD position is 2 and DP position is 3.
  # Different version of GATK tool may assign different names of the variables. Please refer to the GATK manual for more information.
  head(Format)
  GT_pos=1
  AD_pos=2
  DP_pos=3

  # functions to create X, N, and Z vectors
  exinfo_x <- function(info){
    return(sapply(strsplit(info,':'),function(x){ifelse(x[GT_pos]=='./.',NA,as.numeric(strsplit(x[AD_pos],',')[[1]][2]))},simplify=T))
  }
  exinfo_n <- function(info){
    return(sapply(strsplit(info,':'),function(x){ifelse(x[GT_pos]=='./.',NA,as.numeric(x[DP_pos]))},simplify=T))
  }
  exinfo_z <- function(info){
    return(sapply(strsplit(info,':'),function(x){ifelse(x[GT_pos]=='./.',NA,ifelse(x[GT_pos]!='0/0',1,0))},simplify=T))
  }

  # Create X, N, and Z vectors
  ## If these vectors take a long time to create, store dlist and run these lines on HPC cluster instead
  X=sapply(dlist,exinfo_x)
  N=sapply(dlist,exinfo_n)
  Z=sapply(dlist,exinfo_z)

  # Quick filtering
  # If 99% of the cells did not detect any read, we remove such entry.
  thres=floor(0.01*ncol(X))
  sel=rowSums(!is.na(Z))>thres
  X=X[sel,]
  N=N[sel,]
  Z=Z[sel,]
  list(Info=Info, X=X, N=N, Z=Z)
}

#' Create Distance Matrix of SNVs
#'
#' This function creates a distance matrix of SNVs using alternative allele
#' read counts, total allele read counts, and a mutation profile matrix
#'
#' @param filevcf Path to input VCF file
#' @param mutations Path to input mutations RDA file
#' @param SRR_ID Path to file with SRR IDs
#' @param output_csv Name of output csv file with distance matrix
#' @return SNV Distance Matrix
#' @export
createSNVMatrix <- function(filevcf, mutations, SRR_ID, output_csv){
  # Load and parse Mutation Data
  mutations<-vcf_to_SNVinput(filevcf)
  X<-as.matrix(mutations$X);N<-as.matrix(mutations$N);Z<-as.matrix(mutations$Z);
  X[is.na(X)]<-0;N[is.na(N)]<-0;Z[is.na(Z)]<-0;
  cells<-gsub("sampleS", "",unlist(strsplit(colnames(X), "\\."))[1:194*2-1])

  # Extract Cell information and add to Mutation Data
  sampleID<-read.table(SRR_ID, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  sampleID<-sampleID[cells,]; label=sampleID$type
  colnames(X)<-sampleID$Title;colnames(N)<-sampleID$Title;colnames(Z)<-sampleID$Title;
  data<-list(X=X, N=N, Z=Z,label=label)

  # Filter cells and create distance matrix
  data_qc = FilterCellMutation(data$X,data$N,data$Z,data$Info,data$label,cut.off.VAF = 0.05, cut.off.sd = 5)
  data_qc$dist = DENDRO.dist(data_qc$X,data_qc$N,data_qc$Z,show.progress=FALSE)

  # Write distance matrix to csv file
  write.csv(as.matrix(data_qc$dist), file=output_csv)
}
