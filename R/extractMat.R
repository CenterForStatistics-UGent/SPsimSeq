#' A function with S4 dispatching to extract the count matrix
#' @param Y a matrix, data frame or SingleCellExperiment
#' @return A data matrix with samples in the columns and genes in the rows
#' @rdname extractMat
setGeneric("extractMat", function(Y, ...) standardGeneric("extractMat"))
#' @rdname extractMat
#' @import methods
#' @importFrom Biobase exprs
setMethod("extractMat", "SingleCellExperiment", function(Y, ...){
  counts(Y)
})
#' @rdname extractMat
setMethod("extractMat", "matrix", function(Y, ...){
  Y
})
#' @rdname extractMat
setMethod("extractMat", "data.frame", function(Y, ...){
  as.matrix(t(Y))
})