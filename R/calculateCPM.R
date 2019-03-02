#' Calculates counts per millions of reads
#' 
#' @param X an object with class 'SingleCellExperiment', 'data.frame', or 'matrix' 
#' typically containing gene expression data with genes in rows
#' @param const.mult a numerical constant indicatin the desired library size for all samples/cells
#' @param  ... further arguments passed to or from other methods.
#' 
#' @return a matrix of CPM
#' @export
#' @examples  
#' dat <- make.example.data(n.gene = 10,n.sample = 20, n.group = 1, n.batch = 1)
#' cpm.dat <- calCPM(dat)

calCPM <- function(X, const.mult=1e6, ...){
  if(class(X)=="SingleCellExperiment"){
    x <- counts(X)
    if(all(dim(x)>1)){
      as.matrix(x) %*% diag(1/colSums(as.matrix(x))) * const.mult
    }
    else{
      stop("Calculating CPM for a non matrix data!")
    }
  }
  else if(class(X) %in% c("data.frame", "matrix")){
    if(all(dim(X)>1)){
      as.matrix(X) %*% diag(1/colSums(as.matrix(X))) * const.mult
    }
    else{
      stop("Calculating CPM for a non matrix data!")
    }
  }
  else{
    stop("Calculating CPM for a data with class not in 'SingleCellExperiment', 'data.frame', 'matrix'!")
  }  
}