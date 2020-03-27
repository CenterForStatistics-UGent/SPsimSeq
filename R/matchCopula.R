#' Match copulas to estimated SP distribution
#'
#' @param cumDens The cumulative densities evaluated
#' @param exprmt.design the design
#' @param copSam the sampled copula
#' @param sel.genes.ii the gene
#'
#' @return the outcome values as a vector
#' @importFrom stats runif
matchCopula = function(cumDens, exprmt.design, copSam, sel.genes.ii){
  vapply(seq_along(exprmt.design$sub.batchs), function(i){
    copula = copSam[sel.genes.ii,i]
    cd = cumDens[[i]]
    id = which.min(abs(cd$Gy-copula))
    runif(1, cd$breaks[id], cd$breaks[id+1])
  }, FUN.VALUE = numeric(1))
}