#' A function to obtain copulas or uniform random variables
#'
#' @param cpm.data the transformed data matrix
#' @param batch the batch indicators
#' @param n.batch number of batches
#'
#' @return The estimated correlation matrices per batch
obtCorMatsBatch <- function(cpm.data, batch, n.batch){
  lapply(sort(unique(batch)), function(bb){ 
    voom.res <- limma::voom(X[, batch==bb])
    W <- voom.res$weights
    WGCNA::cor(x=t(X[, batch==bb]), weights.x = t(W), 
               use = "pairwise.complete.obs")
  }) 
}

Z       <- mvtnorm::rmvnorm(n=n.batch[bb], sigma = cor.mat)  
Cpl     <- apply(Z, 2, function(z) pnorm(z))
colnames(Cpl) <- rownames(X)
Cpl