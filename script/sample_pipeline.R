# To install CAISC: devtools::install_github("lizamathews/CAISC")
# The following code shows a sample of how to use the commands in the CAISC package

createSNVMatrix("GSE57872_merged.vcf", "DENDRO_input_10mutations.rda", "SRR_ID.txt", "snvDistanceMatrix.csv")
createCNVMatrix("expressionmatrix.txt", "annotationCell.txt", "annotationGene.txt", "GSE57872.RData")

#load
snv_matrix <- read.csv("snvDistanceMatrix.csv")
cnv_matrix <- read.csv("CNVdistMatrix.zscore.csv")

if (identical(dim(snv_matrix), dim(cnv_matrix))){
  integrated_matrix <- createIntegratedMatrix(cnv_matrix, snv_matrix)
}



write.csv(integrated_matrix, file="integratedMatrix.csv", quote = FALSE)

