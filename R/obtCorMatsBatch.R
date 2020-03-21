#' A function to obtain copulas or uniform random variables
#'
#' @param cpm.data the transformed data matrix
#' @param batch the batch indicators
#' @param n.batch number of batches
#'
#' @return The estimated correlation matrices per batch
obtCorMatsBatch <- function(cpm.data, batch, n.batch){
  lapply(sort(unique(batch)), function(bb){ 
    voom.res <- limma::voom(cpm.data[, batch==bb])
    W <- voom.res$weights
    WGCNA::cor(x=t(cpm.data[, batch==bb]), weights.x = t(W), 
               use = "pairwise.complete.obs")
  }) 
}
