# Calculates counts per millions of reads
calCPM <- function(X, const.mult=1e6, norm.lib.size=TRUE, norm.factors=NULL, 
                   logt=FALSE, log_base=2, prior.count=1, ...){
  if(is(X, "SingleCellExperiment") | is(X, "matrix") | is(X, "data.frame")){
    if(!all(dim(X)>1)){
      stop("Calculating CPM for a non matrix data!")
    }
    else if(logt & (log_base%%1 != 0 | log_base<2 | prior.count<0)){
      stop("Invalid log base or prior count!")
    }
    else if(norm.lib.size & !is.null(norm.factors)){
      if(length(norm.factors) != ncol(X) | !is.numeric(norm.factors)){
        stop("Invalid normalization factors passed!")
      }
    }
  }else{
    stop("Calculating CPM for a data with class not in 
             'SingleCellExperiment', 'data.frame', 'matrix'!")
  }
  
  if(is(X, "SingleCellExperiment")){
    x <- counts(X)
  }
  else{
    x <- X
  }
  
  if(norm.lib.size & is.null(norm.factors)){
    nf <- edgeR::calcNormFactors(x)
  }
  else if(norm.lib.size & !is.null(norm.factors)){
    nf <- norm.factors
  }
  else{
    nf <- rep(1, ncol(x))
  }
  cpm <- as.matrix(x) %*% diag(1/(nf*colSums(as.matrix(x)))) * const.mult
  
  if(logt){
    lcpm <- log(cpm+prior.count, base=log_base)
    lcpm
  }
  else{
    cpm
  }
}