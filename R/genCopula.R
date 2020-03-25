#' Generate a copula instance
#' @param corMat list of correlation matrices
#' @param exprmt.design Number of batches, and batch vector
#' 
#' @return a list of copula instances
#' @importFrom mvtnorm rmvnorm
genCopula = function(corMats, exprmt.design){
  n.batch = rowSums(exprmt.design$exprmt.config)
  copList = lapply(names(n.batch), function(b){
    Z       <- rmvnorm(n = n.batch[b], sigma = corMats[[b]])  
    Cpl     <- apply(Z, 2, pnorm)
    colnames(Cpl) <- rownames(corMats[[b]])
    t(Cpl)
  })
  Reduce(copList, f = cbind)
}