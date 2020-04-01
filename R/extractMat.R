#' A function with S4 dispatching to extract the count matrix
#' @param Y a matrix, data frame, phyloseq object or SingleCellExperiment
#' @param ... additional arguments, currently ignored
#' @return A data matrix with samples in the columns and genes in the rows
#' @rdname extractMat
setGeneric("extractMat", function(Y, ...) standardGeneric("extractMat"))
#' @rdname extractMat
#' @import methods
#' @importFrom SingleCellExperiment counts
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
#' @rdname extractMat
#' @importFrom phyloseq t taxa_are_rows otu_table
setMethod("extractMat", "phyloseq", function(Y, ...){
  if(!taxa_are_rows(Y)){
    Y = t(Y)
  }
  as(otu_table(Y), "matrix")
})