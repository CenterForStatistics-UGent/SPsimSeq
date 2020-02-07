# Fit log linear model for each gene
# 
# @param yy a list object contating a result from obtCount() function for a single gene 
# @param  ... further arguments passed to or from other methods.
# 
# @return a list object contating the fitted log linear model
# @examples 
# # 
# @importFrom stats glm pnorm coef vcov
fitLLmodel <- function(yy, ...){
  Ny=yy$Ny
  x=yy$S
  ny=sum(Ny)
  g0= suppressWarnings(pnorm(yy$uls, yy$mu.hat, yy$sig.hat)) - 
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
  if(!is.null(llm)){
    parm.list=list(betas=coef(llm), v=stats:::vcov.glm(llm), 
                   mu.hat=yy$mu.hat, sig.hat=yy$sig.hat)
  }
  else{
    parm.list=list(betas=NULL, v=NULL, mu.hat=yy$mu.hat, sig.hat=yy$sig.hat)
  }
  list(yy=yy, llm=llm, parm.list=parm.list)
}
#' Fast fit Poisson regression
#'
#' @param Ny vector of counts
#' @param x regressor
#' @param degree degree of the polynomial
#' @param offset offset
#'
#' @return see glm.fit
fitPoisGlm = function(Ny, x, degree, offset){
  desMat = buildXmat(x, degree+1)
  glm.fit(y = Ny, x = desMat, family = poisson(), offset = offset)
}
