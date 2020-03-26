#' A function that generates the simulated data for a single gene
#'
#' @param cumDens cumulative density
#' @param sel.genes.ii selected gene
#' @param exprmt.config experiment configuration
#' @param log.CPM.transform a boolean, is log-CPM transform required?
#' @param prior.count the prior count
#' @param LL the library sizes
#' @param copSam the generated copula
#' @param model.zero.prob a boolean, should the zeroes be modelled separately
#' @param fracZero.logit.list The zero model
#' @param const.mult a large constant for the CPM transform, normally 1e6
#'
#' @return Simulated cpm values
SPsimPerGene <- function(cumDens, exprmt.design, 
                         sel.genes.ii, log.CPM.transform, prior.count, LL,
                         copSam, model.zero.prob, fracZero.logit.list, 
                         const.mult){
    ## Match with copula to simulate data 
    Y.star <- matchCopula(cumDens = cumDens, exprmt.design = exprmt.design, 
                          copSam = copSam, sel.genes.ii = sel.genes.ii)
    ## Back tranform to counts if needed
    if(log.CPM.transform){
      Y.star = round(((exp(Y.star) - prior.count)*LL)/const.mult)
      Y.star[Y.star<0] = 0L
    }
    ## Add zeroes if needed
    if(model.zero.prob){
      logLL = log(LL)
      zeroIds = vapply(seq_along(Y.star), FUN.VALUE = logical(1), function(i){
        samZeroID(fracZero.logit.list[[exprmt.design$sub.batchs[[i]]]], 
                   logLL = logLL[i], gene = sel.genes.ii)
      })
      Y.star[zeroIds] = 0L
    }
   return(Y.star)
}
