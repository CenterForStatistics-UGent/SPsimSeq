#' extract the glm parameters
#'
#' @param llm 
#'
#' @return a list of regression parameters and variance covariance matrix
extractGlmParams = function(llm){
  if(!is.null(llm)){
    class(llm) = "glm"
    betas = coef(llm)
    v = vcov(llm)
  }
  else{
    betas = numeric(5)
    v = matrix(0, 5, 5)
  }
  list(betas = betas, v = v)
}