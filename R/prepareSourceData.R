#' Prepare source data for simulation.
#' 
#' @param s.data a source data matrix
#' @param batch a vector containg btach indicator for each sample/cell
#' @param group a vector containg group indicator for each sample/cell
#' @param cand.DE.genes a list object contating candidate predictor (DE) genes
#' @param exprmt.design a list contating simulation design configuration 
#' @param lfc.thrld,llStat.thrld,t.thrld,w see ?chooseCandGenes
#'
#' @return a list object
prepareSourceData <- function(s.data, batch, group, cand.DE.genes,  
                              exprmt.design, lfc.thrld, llStat.thrld,
                              t.thrld, w, log.CPM.transform, prior.count, const.mult){
  # design element
  n.batch <- exprmt.design$n.batch
  n.group <- exprmt.design$n.group
  config.mat <- exprmt.design$exprmt.config
  
  # calculate log CPM 
  cpm.data <- calculateCPM(s.data, log.CPM.transform = log.CPM.transform, 
                          prior.count = prior.count, const.mult = const.mult)
  # subset batches
  if(!is.null(batch)){ 
    if(length(n.batch) < length(unique(batch))){
      sub.batchs <- sort(sample(length(unique(batch)), length(n.batch))) 
    }else if(length(n.batch) == length(unique(batch))){
      sub.batchs <- seq_along(n.batch)
    }
  }else if(is.null(batch)){
    sub.batchs <- 1
    batch <- rep(1, ncol(s.data))
  }
  # select candidate genes
  if(is.null(cand.DE.genes)){
    cand.DE.genes = if(!is.null(group) & length(unique(group))>1){ 
     chooseCandGenes(cpm.data = cpm.data, X = group, const = prior.count, 
                     lfc.thrld = lfc.thrld, t.thrld = t.thrld, 
                     llStat.thrld = llStat.thrld, w = w)
    } else {
      list(null.genes= rownames(s.data)) 
    } 
  }
  list(cand.DE.genes=cand.DE.genes, cpm.data=cpm.data, sub.batchs=sub.batchs)
}
