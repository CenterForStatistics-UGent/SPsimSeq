#' Prepare source data for simulation.
#' 
#' @param s.data a source data (a SingleCellExperiment object)
#' @param batch a vector containg btach indicator for each sample/cell
#' @param group a vector containg group indicator for each sample/cell
#' @param cand.genes a list object contating candidate genes
#' @param exprmt.design a list contating simulation design configuration 
#' @param  ... further arguments passed to or from other methods.
#'
#' @return a list object
#'
#' @examples
#'  # example
#' @export 

prepareSourceData <- function(s.data, batch=NULL, group=NULL, cand.genes=NULL,  
                              exprmt.design, ...){
  
  # design element
  n.batch <- exprmt.design$n.batch
  n.group <- exprmt.design$n.group
  config.mat <- exprmt.design$exprmt.config
  
  # calculate log CPM
  if(class(s.data)=="SingleCellExperiment"){
    cpm.data <- log(calCPM(counts(s.data))+1)
    L <- colSums(counts(s.data))
  }
  else if(class(s.data) %in% c("data.frame", "matrix")){
    cpm.data <- log(calCPM(s.data)+1)
    L <- colSums(s.data)
  }
  
  # subset batches
  if(!is.null(batch)){ 
    if(length(n.batch) < length(unique(batch))){
      sub.batchs <- sort(sample(length(unique(batch)), length(n.batch))) 
    }
    else if(length(n.batch) == length(unique(batch))){
      sub.batchs <- 1:length(n.batch)
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
  LL <- lapply(1:length(sub.batchs), function(b){
    L.b <- L[batch==sub.batchs[b]]
    fit.ln <- fitdist(L.b, distr = "lnorm")$estimate
    L.b.pred <- rlnorm(n.batch[b], fit.ln[["meanlog"]], fit.ln[["sdlog"]])
    
    ## randomly split into the groups
    gr <- rep(1:length(n.group), config.mat[b, ])
    split(L.b.pred, gr)
    #LL.splited <- lapply(sort(unique(group)), function(g) L.b.pred[group==g & batch==b])
    #names(LL.splited) <- paste0("grp_", sort(unique(group)))
    #LL.splited
  }) 
  
  
  # select candidate genes
  if(is.null(cand.genes) & !is.null(group) & length(unique(group))>1){ 
    X <- group
    cand.genes <- chooseCandGenes(s.data=s.data, X=X, ...)
  } 
  else if(is.null(cand.genes) & (is.null(group) | length(unique(group))==1)){
    cand.genes <- list(null.genes= rownames(s.data))
    group <- rep(1, ncol(s.data))
  }
  
  list(cand.genes=cand.genes, LL=LL, cpm.data=cpm.data, sub.batchs=sub.batchs)
}
