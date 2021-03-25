#' Calculate entropy
#'
#' This function calculates entropy using the degree.
#'
#'
#' @param freqs
calculateEntropy <- function(freqs) {
  freqs<-freqs/sum(freqs)
  -sum(freqs * log2(freqs))
}

#' Calculate entropy matrix
#'
#' This function calculates a matrix of entropy values given a matrix.
#'
#'
#' @param matrixWeighted
calculateEntropyMatrix <- function(matrixWeighted) {
  matrixWeighted<-(matrixWeighted-min(matrixWeighted))/(max(matrixWeighted)-min(matrixWeighted))
  degreeNode<-rowSums(matrixWeighted)
  calculateEntropy(degreeNode)
}

#' Calculate entropy matrix for all edges
#'
#' This function calculates entropy values for all edges in the given matrix.
#'
#' @param matrixWeighted input matrix
calculateEntropyMatrixAllEdge <- function(matrixWeighted) {
  matrixWeighted<-(matrixWeighted-min(matrixWeighted))/(max(matrixWeighted)-min(matrixWeighted))
  NodeProbability<-array(matrixWeighted/sum(matrixWeighted))
  NodeProbability<-NodeProbability[NodeProbability>0]  #Dialogo are 0
  calculateEntropy(NodeProbability)
}

#' Calculate alpha values
#'
#' This function calculates alpha values given theta and a list of entropy values.
#'
#' @param entropylist list of entropy values
#' @param theta number of networks to be integrated
EntropytoALPHA <- function(entropylist, theta) {
  ALPHA=1-exp(-1/abs(entropylist)^theta)
  ALPHA/sum(ALPHA)
}

#' Create Integrated Matrix
#'
#' This function creates the integrated CNV and SNV cell-cell distance matrix
#' using an entropy weighted methods.
#'
#' @param cnvDistMatrix cell-cell distance matrix generated with CNV data
#' @param snvDistMatrix cell-cell distance matrix generated with SNV data
#' @return distcombined, integrated matrix
#' @export
createIntegratedMatrix <- function(cnvDistMatrix, snvDistMatrix){
  snvDistMatrix<-(snvDistMatrix-mean(snvDistMatrix))/sd(snvDistMatrix)
  overlapCells<-intersect(row.names(snvDistMatrix), rownames(cnvDistMatrix))
  snvDistMatrix<-snvDistMatrix[overlapCells,overlapCells]
  cnvDistMatrix<-cnvDistMatrix[overlapCells, overlapCells]

  calculateEntropyMatrix(snvDistMatrix)
  calculateEntropyMatrix(cnvDistMatrix)
  entropylist<-c(calculateEntropyMatrix(snvDistMatrix),calculateEntropyMatrix(cnvDistMatrix))
  alphavalues<-EntropytoALPHA(entropylist,2)
  distcombined<-alphavalues[1]*snvDistMatrix + alphavalues[2]*cnvDistMatrix
  return(distcombined)

}
