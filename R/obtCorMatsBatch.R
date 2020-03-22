#' A function to obtain copulas or uniform random variables
#'
#' @param cpm.data the transformed data matrix
#' @param batch the batch indicators
#'
#' @return The estimated correlation matrices per batch
obtCorMatsBatch <- function(cpm.data, batch){
  tapply(colnames(cpm.data), batch,  function(coln){ 
    data = cpm.data[, coln]
    voom.res <- limma::voom(data)
    W <- voom.res$weights
    WGCNA::cor(x=t(data), weights.x = t(W), use = "pairwise.complete.obs")
  }) 
}
