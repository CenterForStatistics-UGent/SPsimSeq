#' Fit log linear model for each gene
#' 
#' @param yy a list object contating a result from obtCount() function for a single gene 
#' @param mu.hat,sig.hat Carrier density estimators
#' @param n number of observations
#' 
#' @return a list object containing the fitted log linear model and carrier density
#' @importFrom stats pnorm
fitLLmodel <- function(yy, mu.hat, sig.hat, n){
  #Evaluate carrier density
  g0 = diff(pnorm(yy$breaks, mean = mu.hat, sd = sig.hat))
  ofs = log(g0*n)
  llm = NULL; degree = 4
  while(is.null(llm) && (degree >= 1)){
  llm <- tryCatch(fitPoisGlm(yy$counts, yy$mids, degree, offset = ofs), 
                  error=function(e){}, warning=function(w){}) 
  if(!is.null(llm) && llm$rank != ncol(llm$R)){
    llm = NULL
  }
  degree = degree - 1
  }
  return(c(yy, list(g0 = g0, coef = llm$coef)))
}
