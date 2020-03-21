#' Calculates counts per millions of reads, possibly with log-transform
#'
#' @param X raw data matrix
#' @param const.mult a constant to multiply with
#' @param log.CPM.transform a boolean, is log-transform desired
#' @param prior.count prior count to be added to the zeroes
#'
#' @return a normalized data matrix
#' @importFrom edgeR calcNormFactors
calulateCPM <- function(X, const.mult, log.CPM.transform, 
                        log_base, prior.count){
  if(log.CPM.transform){
    norm.factors = edgeR::calcNormFactors(X)*colSums(X)
    cpm <- X %*% diag(1/(norm.factors)) * const.mult
    log(cpm+prior.count)
  } else X
}