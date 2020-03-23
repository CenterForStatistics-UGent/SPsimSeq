#' Gene level param estimates for density estimation
#'
#' @param cpm.data.i full vector of genewise observation
#' @param batch,group batch and group information 
#' @param de.ind a boolean, is the gene to be DE?
#' @param model.zero.prob a boolean, should zero-density be modelled?
#' @param min.val minimum value
#' @param w weight
#'
#' @return list of density estimates
geneParmEst <- function(cpm.data.i, batch, group, 
                          de.ind, model.zero.prob, min.val, w){
  tapply(seq_along(cpm.data.i), batch, function(i){ 
    #When not DE, overal density estimation
    if(!de.ind){
      parmEstDensVec(cpm.data.i[i], model.zero.prob, min.val, w = w)
    } else {
      #Otherwise separately per group
     tapply(cpm.data.i[i], group[i], function(Y){
       parmEstDensVec(Y, model.zero.prob, min.val, w = w)
     }) 
    }
  })
}
#' Density estimation on a single vector
#'
#' @param Y0 the vector of observations
#' @param model.zero.prob,min.val,w see geneParmEst()
#' @param prev.min.val minimum prevalence of minimum values
#' @param min.count.nonnull minimum count for estimation
#'
#' @return density estimates
parmEstDensVec = function(Y0, model.zero.prob, min.val, w, prev.min.val = 0.25, 
                           min.count.nonnull = 3){
  #Subset vector
  Y = if(model.zero.prob & mean(Y0==min.val) > prev.min.val){
    Y0[Y0>min.val]
  } else Y0
if(sum(Y > min.val) > min.count.nonnull){
  #Fit normal carrier density
  mu.hat = mean(Y)
  sig.hat = sd(Y)
  #Bin counts
  countY <- obtCount(Y = Y, w = w)
  #Fit exponential density
  fitLLmodel(yy = countY, mu.hat = mu.hat, sig.hat = sig.hat, n = length(Y))
} 
else{
  list(g0 = diff(pnorm(countY$breaks, mean = mu.hat, sd = sig.hat)))
  }
}