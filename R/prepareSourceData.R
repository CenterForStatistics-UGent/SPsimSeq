#' Prepare source data for simulation.
#' 
#' @param s.data a source data (a SingleCellExperiment object)
#' @param batch a vector containg btach indicator for each sample/cell
#' @param group a vector containg group indicator for each sample/cell
#' @param cand.genes a list object contating candidate genes
#' @param exprmt.design a list contating simulation design configuration 
#' @param fc.type a character indicating how DE genes should be calculated ('mean.diff'= difference in the mean log CPM, 'rank.sum'= a U-statistic for the log CPM)
#' @param lfc.thrld a numeric value for the minimum fold change for DE genes (if fc.type='mean.diff')
#' @param t.thrld a numeric value for the minimum t statistic for DE genes (if fc.type='mean.diff')
#' @param U.thrld a numeric value for the minimum U-statistic for DE genes (if fc.type='rank.sum')
#' @param llStat.thrld a numeric value for the minimum squared test statistics from a log-linear model containing X as a covariate to select DE genes
#' @param carrier.dist a character indicating the type of carrier density (carrier.dist="normal" or carrier.dist="kernel")
#' @param w a numeric value between 0 and 1 or NULL refering the number of classes to be created for the outcome data (if NULL the algorithm to calculate breakes in graphics::hist() function will be used)
#' @param max.frac.zero a numeric value between 0 and 1 indicating the maximum fraction of zero counts that a DE gene should have
#' @param max.frac.zeror.diff a numeric value between 0 and 1 indicating the maximum  absolute difference in the fraction of zero counts between the groups for DE genes
#' @param  ... further arguments passed to or from other methods.
#'
#' @return a list object
#'
#' @examples
#'  # example
#' @export 

prepareSourceData <- function(s.data, batch=NULL, group=NULL, cand.genes=NULL,  
                              exprmt.design, fc.type="mean.diff", lfc.thrld=0, U.thrld=0.7,
                              llStat.thrld=5, t.thrld=2.5, carrier.dist="normal",
                              w=0.5, max.frac.zero=0.7, max.frac.zeror.diff=1, ...){
  
  # design element
  n.batch <- exprmt.design$n.batch
  n.group <- exprmt.design$n.group
  config.mat <- exprmt.design$exprmt.config
  
  # calculate log CPM
  if(class(s.data)=="SingleCellExperiment"){
    cpm.data <- log(edgeR::cpm(counts(s.data))+1)
    L <- colSums(counts(s.data))
  }
  else if(class(s.data) %in% c("data.frame", "matrix")){
    cpm.data <- log(edgeR::cpm(s.data))
    L <- colSums(s.data)
  }
  
  # subset batches
  if(length(n.batch) < length(unique(batch))){
    sub.batchs <- sort(sample(length(unique(batch)), length(n.batch))) 
  }
  else if(length(n.batch) == length(unique(batch))){
    sub.batchs <- 1:length(n.batch)
  }
  else{
    stop("Invalid number of batches passed!")
  }
  
  # simulate library size
  LL <- lapply(1:length(sub.batchs), function(b){
    L.b <- L[batch==sub.batchs[b]]
    fit.ln <- fitdistrplus::fitdist(L.b, distr = "lnorm")$estimate
    L.b.pred <- rlnorm(n.batch[b], fit.ln[["meanlog"]], fit.ln[["sdlog"]])
    
    ## randomly split into the groups
    gr <- rep(1:length(n.group), config.mat[b, ])
    split(L.b.pred, gr)
  }) 
  
  
  # select candidate genes
  if(is.null(cand.genes) & !is.null(group) & length(unique(group))==2){ 
    X <- group
    cand.genes <- chooseCandGenes(s.data=s.data, X=X, fc.type=fc.type, lfc.thrld=lfc.thrld,
                                  U.thrld=U.thrld, llStat.thrld=llStat.thrld, t.thrld=t.thrld,
                                  carrier.dist=carrier.dist, w=w, max.frac.zero=max.frac.zero)
  } 
  else if(is.null(cand.genes) & is.null(group)){
    cand.genes <- list(null.genes= rownames(s.data))
  }
  
  list(cand.genes=cand.genes, LL=LL, cpm.data=cpm.data, sub.batchs=sub.batchs)
}
