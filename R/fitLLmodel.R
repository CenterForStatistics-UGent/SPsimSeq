#' Fit log linear model for each gene
#' 
#' @param yy a list object contating a result from obtCount() function for a single gene 
#' @param  ... further arguments passed to or from other methods.
#' 
#' @return a list object contating the fitted log linear model
#' @examples 
#' # 
#' @import stats

fitLLmodel <- function(yy, ...){
  Ny=yy$Ny
  x=yy$S
  ny=sum(Ny)
  g0=pnorm(yy$uls, yy$mu.hat, yy$sig.hat) - pnorm(yy$lls, yy$mu.hat, yy$sig.hat)
  ofs=1 
  llm <- try(glm(Ny~I(x)+I(x^2)+I(x^3)+I(x^4), family = "poisson", offset = log(g0*ny+ofs)), 
             silent = TRUE)
  if(class(llm)[1] != "try-error" & (llm$rank != ncol(llm$R))){
    llm <- try(glm(Ny~I(x)+I(x^2)+I(x^3), family = "poisson", offset = log(g0*ny+ofs)), 
               silent = TRUE)
    if(class(llm)[1] != "try-error" &  (llm$rank != ncol(llm$R))){
      llm <- try(glm(Ny~I(x)+I(x^2), family = "poisson", offset = log(g0*ny+ofs)), silent = TRUE)
      if(class(llm)[1] != "try-error" &  llm$rank != ncol(llm$R)){
        llm <- try(glm(Ny~I(x), family = "poisson", offset = log(g0*ny+ofs)), silent = TRUE)
      }
      else{
        llm <- NULL
      }
    } 
  } 
  
  if(!is.null(llm)){
    parm.list=list(betas=coef(llm), v=vcov(llm), mu.hat=yy$mu.hat, sig.hat=yy$sig.hat)
  }
  else{
    parm.list=list(betas=NULL, v=NULL, mu.hat=yy$mu.hat, sig.hat=yy$sig.hat)
  }
  
  list(yy=yy, llm=llm, parm.list=parm.list)
}
