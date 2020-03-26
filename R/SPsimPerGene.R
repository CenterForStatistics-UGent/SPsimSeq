#' A function that generates the simulated data for each gene
#'
#' @param densList.ii 
#' @param DE.ind.ii 
#' @param sel.genes.ii 
#' @param exprmt.config
#' @param log.CPM.transform 
#' @param prior.count
#' @param LL 
#' @param copSam 
#' @param model.zero.prob 
#' @param fracZero.logit.list 
#'
#' @return Simulated cpm values
SPsimPerGene <- function(densList.ii, DE.ind.ii, exprmt.design, 
                         sel.genes.ii, log.CPM.transform, prior.count, LL,
                         copSam, model.zero.prob, fracZero.logit.list, 
                         const.mult){
    ## construct the density
    cumDens <- constructDens(densList.ii = densList.ii, 
                             exprmt.design = exprmt.design, 
                             DE.ind.ii = DE.ind.ii)
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
