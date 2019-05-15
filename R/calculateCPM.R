#' Calculates counts per millions of reads
#' 
#' @param X an object with class 'SingleCellExperiment', 'data.frame', or 'matrix' 
#' typically containing gene expression data with genes in rows
#' @param const.mult a numerical constant indicatin the desired library size for 
#' all samples/cells
#' @param norm.lib.size logical value. If TRUE, normalized library size will be used to
#' calculate CPM. In particulr, TMM normalization factors will be used unless @param 
#' norm.factors is not NULL
#' @param norm.factors a numerical vector of normalization factors
#' @param logt a logical value. If TRUE, log(base=@param log_base) of the CPM 
#' will be returned
#' @param log_base an integer for the base of logarithmic transformation.
#' It will only be considered if @param log = TRUE.
#' @param prior.count a positive integer to be added to the CPM prior to log transformation.
#' It will only be considered if @param log = TRUE.
#' @param  ... further arguments passed to or from other methods.
#' 
#' @return a matrix of CPM
#' @export
#' @examples
#' dat <- make.example.data(n.gene = 10,n.sample = 5, n.group = 1, n.batch = 1)
#' cpm.dat <- calCPM(dat)
#' cpm.dat <- calCPM(dat, norm.lib.size = FALSE)
#' 
#' cpm.dat <- calCPM(dat, logt = TRUE)
#' cpm.dat <- calCPM(dat, logt = TRUE, log_base = 10)
#' 
#' cpm.dat <- calCPM(dat, logt = TRUE)
#' cpm.dat <- calCPM(dat, logt = TRUE, prior.count = 0.25)
#' 
#' @importFrom edgeR calcNormFactors

calCPM <- function(X, const.mult=1e6, norm.lib.size=TRUE, norm.factors=NULL, logt=FALSE,
                   log_base=2, prior.count=1, ...){
  if(class(X) %in% c("SingleCellExperiment", "data.frame", "matrix")){ 
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
  }
  else{
    stop("Calculating CPM for a data with class not in 'SingleCellExperiment', 'data.frame', 'matrix'!")
  }
  
  if(class(X)=="SingleCellExperiment"){
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


