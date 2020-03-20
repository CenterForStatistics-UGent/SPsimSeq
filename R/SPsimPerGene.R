#' A function that generates the simulated data for each gene
#'
#' @param cpm.data 
#' @param est.list.ii 
#' @param DE.ind.ii 
#' @param n.batch 
#' @param n.group 
#' @param batch 
#' @param group 
#' @param config.mat 
#' @param sel.genes.ii 
#' @param log.CPM.transform 
#' @param const 
#' @param min.val 
#' @param null.group 
#' @param LL 
#' @param copulas.batch 
#' @param tot.samples 
#' @param model.zero.prob 
#' @param fracZero.logit.list 
#'
#' @return
SPsimPerGene <- function(cpm.data, est.list.ii, DE.ind.ii, n.batch, n.group, batch, group, config.mat, 
                         sel.genes.ii, log.CPM.transform, const, min.val, null.group, LL,
                         copulas.batch, tot.samples, model.zero.prob, fracZero.logit.list){
    # get batch specific parameters
    par.sample <- obtParSample(est.list.i = est.list.ii, DE.ind.ii = DE.ind.ii, 
                               n.batch = n.batch, group = group)
    # estimate carrier density (g0)
    g0 <- estCarrierDens(est.list.i = est.list.ii, par.sample = par.sample,
                         DE.ind.ii = DE.ind.ii,  group = group)
    # estimate the density g1(y)
    g1 <- estSPDens(est.list.i = est.list.ii, par.sample = par.sample,
                    DE.ind.ii = DE.ind.ii, group = group, g0 = g0)
    # simulate data 
    Y.star <- sampleDatSPDens(cpm.data=cpm.data, sel.genes.i=sel.genes.ii, par.sample=par.sample, 
                              DE.ind.ii=DE.ind.ii, null.group=null.group, LL=LL,
                              copulas.batch=copulas.batch, group=group, batch=batch, 
                              g1=g1, log.CPM.transform=log.CPM.transform, const=const,
                              min.val=min.val, n.group=n.group, n.batch=n.batch,
                              config.mat=config.mat, model.zero.prob=model.zero.prob, 
                              fracZero.logit.list=fracZero.logit.list)
    names(g1) <- names(g0) <- names(par.sample) <- 
      paste0("group_", if(DE.ind.ii) unique(group) else null.group)
     return(Y.star)
}
