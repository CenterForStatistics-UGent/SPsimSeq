# Fit log linear model for each gene
# 
# @param yy a list object contating a result from obtCount() function for a single gene 
# @param  ... further arguments passed to or from other methods.
# 
# @return a list object contating the fitted log linear model
# @importFrom stats glm pnorm coef vcov
fitLLmodel <- function(yy, ...){
  Ny = yy$Ny
  x = yy$S
  ny = sum(Ny)
  g0 = suppressWarnings(pnorm(yy$uls, yy$mu.hat, yy$sig.hat)) - 
    suppressWarnings(pnorm(yy$lls, yy$mu.hat, yy$sig.hat))
  ofs = log(g0*ny+1)
  llm = NULL
  degree = 4
  while(is.null(llm) && (degree >= 1)){
  llm <- tryCatch(fitPoisGlm(Ny, x, degree, offset = ofs), 
                  error=function(e){}, warning=function(w){}) 
  if(!is.null(llm) && llm$rank != ncol(llm$R)){
    llm = NULL
  }
  degree = degree - 1
  }
  est.parms <- parmEstOut(llm)
  return(c(yy, list(betas = est.parms$beta.hat.vec, v = est.parms$V.hat.mat)))
}
