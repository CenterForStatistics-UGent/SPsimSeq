#' Calculates counts per millions of reads, possibly with log-transform
#'
#' @param X raw data matrix
#' @param const.mult a constant to multiply with
#' @param prior.count prior count to be added to the zeroes
#'
#' @return a normalized data matrix
#' @importFrom edgeR calcNormFactors
calculateCPM <- function(X, const.mult, prior.count){
    norm.factors = edgeR::calcNormFactors(X)*colSums(X)
    cpm = X %*% diag(1/(norm.factors)) * const.mult
    colnames(cpm) = colnames(X)
    log(cpm + prior.count)
}