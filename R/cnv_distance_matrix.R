

#library(BiocManager)
#BiocManager::install("infercnv")
#library(infercnv)


#' Create InferCNV Object
#'
#' This function creates an inferCNV object.
#'
#'
#' @param rawCountsMatrix Name of raw counts matrix file
#' @param annotations Name of annotations file
#' @param geneOrder Name of gene order file
#' @param outputFile Name of output file that inferCNV object is saved to
makeInferCNVObject<-function(rawCountsMatrix, annotations, geneOrder, outputFile ){
  infercnv_obj = infercnv::CreateInfercnvObject(raw_counts_matrix= rawCountsMatrix,
                                      annotations_file= annotations,
                                      delim="\t",
                                      gene_order_file= geneOrder,
                                      ref_group_names=c("CSC"))
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir="./",
                               cluster_by_groups=TRUE,
                               analysis_mode = "subclusters",
                               k_obs_groups = 4,
                               denoise=TRUE,
                               HMM=TRUE,
                               HMM_type = "i3",
                               num_threads=6,
                               no_plot=FALSE
  )

  save(infercnv_obj, file=outputFile)
}

#' Create CNV Distance Matrix
#'
#' This function creates a cell-cell distance matrix using CNV data.
#'
#'
#' @param rawCountsMatrix Name of raw counts matrix file
#' @param annotations Name of annotations file
#' @param geneOrder Name of gene order file
#' @param outputFile Name of output file that inferCNV object is saved to
#' @export
createCNVMatrix<-function(rawCountsMatrix, annotations, geneOrder, ref_sample, outputFile ){
  makeInferCNVObject(rawCountsMatrix, annotations, geneOrder, outputFile)
  load(outputFile)
  exp.data<-as.data.frame(infercnv_obj@expr.data)

  # subset for all patient samples, leaving out the reference
  samples<-names(infercnv_obj@tumor_subclusters$subclusters)
  samples<-samples[-(grep(ref_sample, samples, value=FALSE))]

  ###Let us read the cell order
  cellsListDraw<-c()
  for(i in 1:length(samples)){
    cellsList<-infercnv_obj@tumor_subclusters$subclusters[[i]]
    for(ii in 1:length(cellsList)){
      cellsListDraw<-c(cellsListDraw, names(cellsList[[ii]]))
    }
    cat(length(cellsListDraw));cat("\n")
  }



  exp.data.scale<-as.data.frame(t(scale(t(exp.data[,cellsListDraw]))))
  exp.data.scale<-exp.data.scale[!is.na(apply(exp.data.scale, 1,sd)),]

  corValue<-cor(exp.data.scale[, cellsListDraw])
  distValue<-1-corValue
  distValue.zscore<-(distValue-mean(distValue))/sd(distValue)
  distValue.pvalue<-round(1-pnorm(distValue.zscore, lower.tail=TRUE), digits = 6)
  write.csv(corValue, file="corValue.csv", quote = FALSE)
  write.csv(distValue.zscore, file="CNVdistMatrix.zscore.csv", quote = FALSE)
  write.csv(distValue.pvalue, file="CNVdistMatrix.pvalue.csv", quote = FALSE)

}
