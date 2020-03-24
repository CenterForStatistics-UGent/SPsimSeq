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
SPsimPerGene <- function(densList.ii, DE.ind.ii, exprmt.config, 
                         sel.genes.ii, log.CPM.transform, prior.count, LL,
                         copSam, model.zero.prob, fracZero.logit.list, 
                         const.mult){
    ## construct the density
    cumDens <- constructDens(densList.ii = densList.ii, exprmt.config = exprmt.config,
                    DE.ind.ii = DE.ind.ii)
    ## Match with copula to simulate data 
    Y.star <- matchCopula(cumDens = cumDens, exprmt.config = exprmt.config, 
                          copSam = copSam, DE.ind.ii = DE.ind.ii, 
                          sel.genes.ii = sel.genes.ii)
    ## Back tranform to counts if needed
    if(log.CPM.transform){
      Y.star = lapply(seq_along(Y.star), function(i){
        round(((exp(Y.star[[i]]) - prior.count)*LL[[i]])/const.mult)
      })
    }
    ## Add zeroes if needed
    if(model.zero.prob){
      Y.star = lapply(seq_along(Y.star), function(i){
        addZeroes(Y.star[[i]], fracZero.logit.list[[i]])
      })
    }
   return(Y.star)
}
