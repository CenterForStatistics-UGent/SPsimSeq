#' Gene level param estimates for density estimation
#'
#' @param cpm.data.i full vector of genewise observation
#' @param batch,group batch and group information 
#' @param null.group the null group
#' @param sub.batchs 
#' @param de.ind a boolean, is the gene to be DE?
#' @param model.zero.prob a boolean, should zero-density be modelled?
#' @param min.val minimum value
#' @param w weight
#'
#' @return list of density estimates
geneParmEst <- function(cpm.data.i, batch, group, null.group, sub.batchs, 
                          de.ind, model.zero.prob, min.val, w){
  lapply(sub.batchs, function(b){ 
    if(!de.ind){
      Y0 <- cpm.data.i[(batch==b & group==null.group)] 
      parmEstDensVec(Y0, modelzero.prob, min.val)
    } else {
      Y0 <- split(cpm.data.i[batch==b], group[batch==b])
      Y <- lapply(Y0, function(y0){
        parmEstDensVec(y0, modelzero.prob, min.val)
      }) 
    }
  }) 
}
#' Density estimation on a single vector
#'
#' @param Y0 the vector of observations
#' @param model.zero.prob,min.val see geneParmEst()
#' @param prev.min.val minimum prevalence of minimum values
#' @param min.count.nonnull minimum count for estimation
#'
#' @return density estimates
parmEstDensVec = function(Y0, model.zero.prob, min.val, prev.min.val = 0.25, 
                           min.count.nonnull = 3){
  Y = if(model.zero.prob & mean(Y0==min.val) > prev.min.val){
  Y0[Y0>min.val]
  } else Y0
if(sum(Y > min.val) > min.count.nonnull){
  countY   <- obtCount(Y = Y, w = w)
  fitLLmodel(yy = countY)
} 
else{
  list(S=0, lls=0, uls=0, Ny=0, Y=Y, mu.hat=mean(Y), sig.hat=sd(Y), 
       w=w, betas=parmEstOut(llm = NULL)$beta.hat.vec, 
       v = parmEstOut(llm = NULL)$V.hat.mat)
  }
}