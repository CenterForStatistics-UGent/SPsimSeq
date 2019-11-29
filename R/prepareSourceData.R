# Prepare source data for simulation.
# 
# @param s.data a source data (a SingleCellExperiment object)
# @param batch a vector containg btach indicator for each sample/cell
# @param group a vector containg group indicator for each sample/cell
# @param cand.genes a list object contating candidate genes
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

prepareSourceData <- function(s.data, batch=NULL, group=NULL, cand.genes=NULL,  
                              exprmt.design, simCtr, lfc.thrld, llStat.thrld,
                              const, w, log.CPM.transform, lib.size.params=NULL, ...){
  
  # design element
  n.batch <- exprmt.design$n.batch
  n.group <- exprmt.design$n.group
  config.mat <- exprmt.design$exprmt.config
  
  # calculate log CPM 
  if(log.CPM.transform){
    if(is(s.data, "SingleCellExperiment")){ 
      cpm.data <- log(calCPM(counts(s.data))+const)
      L <- colSums(counts(s.data))
    }
    else if(is(s.data, "data.frame") | is(s.data, "matrix")){
      cpm.data <- log(calCPM(s.data)+const)
      L <- colSums(s.data)
    }
  }
  else{
    if(is(s.data, "SingleCellExperiment")){ 
      cpm.data <- counts(s.data)
      L <- colSums(counts(s.data))
    }
    else if(is(s.data, "data.frame") | is(s.data, "matrix")){
      cpm.data <- s.data
      L <- colSums(s.data)
    }
  }
  
  
  # subset batches
  if(!is.null(batch)){ 
    if(length(n.batch) < length(unique(batch))){
      set.seed(simCtr)
      sub.batchs <- sort(sample(length(unique(batch)), length(n.batch))) 
    }
    else if(length(n.batch) == length(unique(batch))){
      sub.batchs <- seq_len(length(n.batch))
    }
    else{
      stop("Invalid number of batches passed: length(n.batch) > length(unique(batch))")
    }
  } 
  else if(is.null(batch)){
    sub.batchs <- 1
    batch <- rep(1, ncol(s.data))
  }
  else{
    stop("Invalid number of batches passed!")
  }
  
  # simulate library size 
  LL <- lapply(seq_len(length(sub.batchs)), function(b){
    if(is.null(lib.size.params)){
      L.b      <- L[batch==sub.batchs[b]]
      fit.ln   <- fitdist(as.numeric(L.b), distr = "lnorm")$estimate 
      L.b.pred <- rlnorm(n.batch[b], fit.ln[["meanlog"]], fit.ln[["sdlog"]])
    }
    else{
      if(length(lib.size.params) != 2 | is.null(names(lib.size.params))){
        stop("The log-normal parameters for the distribution of library sizes must be submitted in a named vector of size 2. 
             Example, lib.size.params = c(meanlog=10, sdlog=0.2). See also ?rlnorm()")
      }else{
        L.b.pred <- rlnorm(n.batch[b], lib.size.params[["meanlog"]], lib.size.params[["sdlog"]])
      } 
    }
    
    
    ## randomly split into the groups
    gr <- rep(seq_len(length(n.group)), config.mat[b, ])
    split(L.b.pred, gr) 
  }) 
  
  
  # select candidate genes
  if(is.null(cand.genes) & !is.null(group) & length(unique(group))>1){ 
    X <- group
    cand.genes <- chooseCandGenes(cpm.data=cpm.data, X=X, const=const,
                                  lfc.thrld=lfc.thrld, llStat.thrld=llStat.thrld, w=w, ...)
  } 
  else if(is.null(cand.genes) & (is.null(group) | length(unique(group))==1)){
    cand.genes <- list(null.genes= rownames(s.data))
    group <- rep(1, ncol(s.data))
  }
  
  list(cand.genes=cand.genes, LL=LL, cpm.data=cpm.data, sub.batchs=sub.batchs)
}
