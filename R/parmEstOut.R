parmEstOut <- function(llm){
  beta.hat.vec <- numeric(5)
  V.hat.mat <- matrix(0, 5, 5)
  
  if(!is.null(llm)){
    class(llm) = "glm"
    betas=coef(llm)
    v=vcov(llm)
  }
  else{
    betas= numeric(5)
    v=matrix(0, 5, 5)
  }
  beta.hat.vec[seq_along(betas)] <- betas
  V.hat.mat[seq_len(nrow(v)), seq_len(ncol(v))] <- v
  
  list(beta.hat.vec=beta.hat.vec, V.hat.mat=V.hat.mat)
}