#' A function to obtain copulas or uniform random variables
#'
#' @param cpm.data the transformed data matrix
#' @param batch the batch indicators
#' @param genewiseCor logical, whther gene wise correlation needs to be simulated
#'
#' @return The estimated correlation matrices per batch
#' @importFrom WGCNA cor
#' @importFrom limma voom
obtCorMatsBatch <- function(cpm.data, batch, genewiseCor){
  tapply(colnames(cpm.data), batch,  function(coln){ 
    if(genewiseCor){
      data <- cpm.data[, coln]
      voom.res <- limma::voom(data)
      W <- voom.res$weights
      WGCNA::cor(x=t(data), weights.x = t(W), use = "pairwise.complete.obs")
    }else{
      data <- cpm.data[, coln]
      cormat <- diag(x=1, nrow = nrow(data))
      rownames(cormat) <- colnames(cormat) <- rownames(data)
      cormat
    } 
  }) 
}
