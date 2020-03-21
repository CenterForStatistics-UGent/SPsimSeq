#' Generate a copula instance
#' @param corMat list of correlation matrices
#' @param n.batch,batch Number of batches, and batch vector
#' 
#' @return a list of copula instances
genCopula = function(corMats, n.batch, batch){
  lapply(sort(unique(batch)), function(bb){ 
    Z       <- mvtnorm::rmvnorm(n=n.batch[bb], sigma = corMats[[bb]])  
    Cpl     <- apply(Z, 2, pnorm)
    colnames(Cpl) <- rownames(corMats[[bb]])
    Cpl
  })
}