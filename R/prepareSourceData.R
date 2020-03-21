#' Prepare source data for simulation.
#' 
#' @param s.data a source data matrix
#' @param batch a vector containg btach indicator for each sample/cell
#' @param group a vector containg group indicator for each sample/cell
#' @param cand.DE.genes a list object contating candidate predictor (DE) genes
#' @param exprmt.design a list contating simulation design configuration 
#' @param const a small constant (>0) to be added to the CPM before log transformation, to avoid  log(0)
#' @param lfc.thrld,llStat.thrld,t.thrld,w see ?chooseCandGenes
#'
#' @return a list object
#'
#' @examples
#'  # example 
prepareSourceData <- function(s.data, batch, group, cand.DE.genes,  
                              exprmt.design, const, lfc.thrld, llStat.thrld,
                              t.thrld, w, log.CPM.transform){
  
  # design element
  n.batch <- exprmt.design$n.batch
  n.group <- exprmt.design$n.group
  config.mat <- exprmt.design$exprmt.config
  
  # calculate log CPM 
  if(log.CPM.transform){
      cpm.data <- log(calCPM(s.data)+const)
  }
  else{
      cpm.data <- s.data
      LS <- colSums(s.data)
  }
  
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
     chooseCandGenes(cpm.data = cpm.data, X = group, const = const, 
                     lfc.thrld = lfc.thrld, t.thrld = t.thrld, 
                     llStat.thrld = llStat.thrld, w = w)
    } else {
      list(null.genes= rownames(s.data)) 
    } 
  }
  list(cand.DE.genes=cand.DE.genes, cpm.data=cpm.data, sub.batchs=sub.batchs)
}
