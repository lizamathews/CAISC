#' Create SNV Distance Matrix input from VCF File
#'
#' This function finds chromosome information, alternative allele read counts (X),
#'  total read counts (N), and mutation profile matrix (Z) from VCF file of mutations
#'
#' @param cnvDistMatrix cell-cell distance matrix generated with CNV data
#' @param snvDistMatrix cell-cell distance matrix generated with SNV data
#' @param integratedDistMatrix integrated distance matrix with CNV and SNV data
#' @param infercnvObj RDA File with inferCNV object
#' @param geneOrder Name of geneOrder file, file with gene annotations
#' @return interactive complex heatmap
#' @export
createComplexHeatmap <- function(snvDistMatrix, cnvDistMatrix, integratedDistMatrix, infercnvObj, geneOrder){
  # clustering function
  # BiocManager::install("InteractiveComplexHeatmap")
  DENDRO.cluster = function(dist,method='ward.D',plot=TRUE,label=NULL,type="phylogram",...){
    clust=hclust(dist,method=method)
    return(clust)
  }
  # load snv matrix
  snvDistMatrix<-(snvDistMatrix-min(snvDistMatrix))/(max(snvDistMatrix)-min(snvDistMatrix))

  # load cnv matrix
  cnvDistMatrix<-(cnvDistMatrix-min(cnvDistMatrix))/(max(cnvDistMatrix)-min(cnvDistMatrix))
  samples<-gsub("\\.", "_", gsub(".bam", "",rownames(cnvDistMatrix)))
  rownames(cnvDistMatrix)<-samples;colnames(cnvDistMatrix)<-samples;

  # find cells that overlap between snv and cnv matrices
  overlapCells<-intersect(row.names(snvDistMatrix), rownames(cnvDistMatrix))
  snvDistMatrix<-snvDistMatrix[overlapCells,overlapCells]
  cnvDistMatrix<-cnvDistMatrix[overlapCells,overlapCells]

  # cluster CNV
  clusters_CNV = DENDRO.cluster(as.dist(cnvDistMatrix),label=rownames(cnvDistMatrix),type='fan', plot = FALSE)
  membershipresult_CNV<-dendextend::cutree(as.dendrogram(clusters_CNV), k = 4)
  RowSideColors_CNV<-rep("black", length(membershipresult_CNV))
  RowSideColors_CNV[membershipresult_CNV==1]<-"red"
  RowSideColors_CNV[membershipresult_CNV==2]<-"brown"
  RowSideColors_CNV[membershipresult_CNV==3]<-"yellow"

  # cluster SNV
  clusters_SNV = DENDRO.cluster(as.dist(snvDistMatrix),label=rownames(cnvDistMatrix),type='fan', plot = FALSE)
  clusters_SNV$labels[clusters_SNV$order]
  membershipresult_SNV<-dendextend::cutree(as.dendrogram(clusters_SNV), k = 4)
  RowSideColors<-rep("black", length(membershipresult_SNV))
  RowSideColors[membershipresult_SNV==1]<-"red"
  RowSideColors[membershipresult_SNV==2]<-"brown"
  RowSideColors[membershipresult_SNV==3]<-"yellow"

  # load cnv object
  load(infercnvObj)
  exp.data<-as.data.frame(infercnv_obj@expr.data)
  samples<-gsub("\\.", "_", gsub(".bam", "",colnames(exp.data)))
  colnames(exp.data)<-samples;

  # load gene order data
  geneLocation<-read.table(geneOrder, row.names = 1)
  geneLocation<-geneLocation[gtools::mixedorder(geneLocation$V2),]
  exp.data<-na.omit(exp.data[row.names(geneLocation), ])
  chroms<-geneLocation[rownames(exp.data),"V2"]
  ColSideColor<-rep("black", dim(exp.data)[1])
  ColSideColor[chroms %in% paste0("chr", 1:11*2)]<-"grey"

  # clustering
  ####Two groups of P2
  cellsListDraw<-clusters_SNV$labels[clusters_SNV$order]   ###order on CNV
  RowSideColors<-RowSideColors[clusters_SNV$order]  #order by clustering

  ####Two groups of P2
  RowSideColors_CNV<-RowSideColors_CNV[clusters_SNV$order]  #order by clustering

  ####By sampleID
  RowSideColors_sample<-rep("black", length(membershipresult_SNV))
  RowSideColors_sample[grepl("P2_1", names(membershipresult_SNV))]<-"red"
  RowSideColors_sample[grepl("P2_2", names(membershipresult_SNV))]<-"yellow"
  RowSideColors_sample<-RowSideColors_sample[clusters_SNV$order]  #order by clustering

  ####clustering result
  dataClustering<-data.frame(RowSideColors_sample=RowSideColors_sample, RowSideColors_CNV=RowSideColors_CNV, RowSideColors_SNV=RowSideColors)
  rownames(dataClustering)<-cellsListDraw

  # specify subset of heatmap data
  subset<-floor(dim(exp.data)[1]/5)

  # create annotation column for heatmap plot
  annotation_column<-as.data.frame(as.matrix(ColSideColor[1:subset*5]))
  rownames(annotation_column)<-rownames(exp.data)[1:subset*5]

  # create data for heatmap plot
  dataPlot<-t(as.matrix(exp.data[1:subset*5,cellsListDraw]))

  # load integrated matrix and make integrated annotations
  clusters_BOTH = DENDRO.cluster(as.dist(integratedDistMatrix),label=rownames(cnvDistMatrix),type='fan', plot = FALSE)
  membershipresult_BOTH<-dendextend::cutree(as.dendrogram(clusters_BOTH), k = 4)
  RowSideColors_BOTH<-rep("black", length(membershipresult_BOTH))
  RowSideColors_BOTH[membershipresult_BOTH==1]<-"red"
  RowSideColors_BOTH[membershipresult_BOTH==2]<-"brown"
  RowSideColors_BOTH[membershipresult_BOTH==3]<-"yellow"

  names(RowSideColors_BOTH)<-names(membershipresult_BOTH)
  dataClustering$RowSideColors_BOTH<-RowSideColors_BOTH[rownames(dataClustering)]
  colnames(dataClustering)<-gsub("RowSideColors_", "", colnames(dataClustering))

  # define annotation colors
  annotation_colors_mergeCNV = list(SNV=c(black="blue",red="red",brown="purple3",yellow="lavenderblush2")
                                    ,CNV=c(black="blue",red="red",brown="purple3",yellow="lavenderblush2")
                                    ,BOTH=c(black="blue",red="red",brown="purple3",yellow="lavenderblush2")
                                    ,sample=c(black="blue",red="red",brown="purple3",yellow="lavenderblush2")
                                    ,V1=c(black="mediumpurple", grey="thistle1"))

  # create interactive complex heatmap
  ht1 = ComplexHeatmap::pheatmap(dataPlot, col = colorRampPalette(c("blue","oldlace","red", "red"))(85), show_rownames	= TRUE, show_colnames = FALSE, legend = FALSE, annotation_legend = FALSE
                                 , cluster_rows=FALSE, cluster_cols=FALSE, annotation_row = dataClustering[, c(3,2)], annotation_col = annotation_column, annotation_colors = annotation_colors_mergeCNV)

  ht1

}
