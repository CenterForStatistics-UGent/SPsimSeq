#' Generate a copula instance
#'
#' @param exprmt.design Number of batches, and batch vector
#' @param corMats List of correlation matrices
#'
#' @return a list of copula instances
#' @importFrom mvtnorm rmvnorm
genCopula = function(corMats, exprmt.design){
  n.batch = rowSums(exprmt.design$exprmt.config)
  copList = lapply(names(n.batch), function(b){
    Z       <- rmvnorm(n = n.batch[b], sigma = corMats[[b]])
    Cpl     <- pnorm(Z)
    colnames(Cpl) <- rownames(corMats[[b]])
    t(Cpl)
  })
  Reduce(copList, f = cbind)
}