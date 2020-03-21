#' Fit log linear model for each gene
#' 
#' @param yy a list object contating a result from obtCount() function for a single gene 
#' @param  ... further arguments passed to or from other methods.
#' 
#' @return a list object containing the fitted log linear model and carrier density
#' @importFrom stats pnorm
fitLLmodel <- function(yy, ...){
  #Carrier density
  g0 = suppressWarnings(pnorm(yy$uls, yy$mu.hat, yy$sig.hat) - 
    pnorm(yy$lls, yy$mu.hat, yy$sig.hat))
  ofs = log(g0*sum(yy$Ny)+1)
  llm = NULL
  degree = 4
  while(is.null(llm) && (degree >= 1)){
  llm <- tryCatch(fitPoisGlm(yy$Ny, yy$S, degree, offset = ofs), 
                  error=function(e){}, warning=function(w){}) 
  if(!is.null(llm) && llm$rank != ncol(llm$R)){
    llm = NULL
  }
  degree = degree - 1
  }
  est.parms <- extractGlmParams(llm)
  return(c(yy, est.parms, list(g0 = g0)))
}
