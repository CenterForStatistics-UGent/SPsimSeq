#' An auxialiary function to quickly construct the polyomial matrix, 
#' using Horner's rule
#'
#' @param x The base 
#' @param nc the number of columns
#'
#' @return A matrix with increasing powers of x in the columns
buildXmat = function(x, nc){
  x.mat <- matrix(1, nrow = length(x), ncol = nc)
  for(k in 2:nc){
    x.mat[,k:nc] <- x.mat[,k:nc]*x
  }
  x.mat
}