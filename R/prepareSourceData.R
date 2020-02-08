# Prepare source data for simulation.
# 
# @param s.data a source data (a SingleCellExperiment object)
# @param batch a vector containg btach indicator for each sample/cell
# @param group a vector containg group indicator for each sample/cell
# @param cand.DE.genes a list object contating candidate predictor (DE) genes
# @param exprmt.design a list contating simulation design configuration 
# @param const a small constant (>0) to be added to the CPM before log transformation, to avoid  log(0).
# default is 1e-5
# @param simCtr (integer) seed number
# @param  ... further arguments passed to or from other methods.
#
# @return a list object
#
# @examples
#  # example 
# @importFrom fitdistrplus fitdist
prepareSourceData <- function(s.data, batch, group, cand.DE.genes,  
                              exprmt.design, simCtr, lfc.thrld, llStat.thrld,
                              const, w, log.CPM.transform, ...){
  
  # design element
  n.batch <- exprmt.design$n.batch
  n.group <- exprmt.design$n.group
  config.mat <- exprmt.design$exprmt.config
  
  # calculate log CPM 
  if(log.CPM.transform){
    if(is(s.data, "SingleCellExperiment")){ 
      cpm.data <- log(calCPM(counts(s.data))+const)
      #LS <- colSums(counts(s.data))
    }else if(is(s.data, "data.frame") | is(s.data, "matrix")){
      cpm.data <- log(calCPM(s.data)+const)
      #LS <- colSums(s.data)
    }
  }
  else{
    if(is(s.data, "SingleCellExperiment")){ 
      cpm.data <- counts(s.data)
      LS <- colSums(counts(s.data))
    }else if(is(s.data, "data.frame") | is(s.data, "matrix")){
      cpm.data <- s.data
      LS <- colSums(s.data)
    }
  }
  
  
  # subset batches
  if(!is.null(batch)){ 
    if(length(n.batch) < length(unique(batch))){
      set.seed(simCtr)
      sub.batchs <- sort(sample(length(unique(batch)), length(n.batch))) 
    }else if(length(n.batch) == length(unique(batch))){
      sub.batchs <- seq_len(length(n.batch))
    }else{
      stop("Invalid number of batches passed: length(n.batch) > length(unique(batch))")
    }
  }else if(is.null(batch)){
    sub.batchs <- 1
    batch <- rep(1, ncol(s.data))
  }else{
    stop("Invalid number of batches passed!")
  }
  
  # # simulate library size  
  # LL <- simLibSize(sub.batchs=sub.batchs, L=LS, lib.size.params = lib.size.params, 
  #                  n.batch=n.batch, batch = batch, n.group = n.group, 
  #                  config.mat = config.mat, ...)
  
  # select candidate genes
  if(is.null(cand.DE.genes) & !is.null(group) & length(unique(group))>1){ 
    X <- group
    cand.DE.genes <- chooseCandGenes(cpm.data=cpm.data, X=X, const=const,
                                  lfc.thrld=lfc.thrld, llStat.thrld=llStat.thrld, w=w, ...)
  }else if(is.null(cand.DE.genes) & (is.null(group) | length(unique(group))==1)){
    cand.DE.genes <- list(null.genes= rownames(s.data)) 
  } 
  list(cand.DE.genes=cand.DE.genes, cpm.data=cpm.data, sub.batchs=sub.batchs)
  #list(cand.DE.genes=cand.DE.genes, LL=LL, cpm.data=cpm.data, sub.batchs=sub.batchs)
}
