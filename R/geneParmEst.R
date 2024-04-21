#' Gene level param estimates for density estimation
#'
#' @param cpm.data.i full vector of genewise observation
#' @param de.ind a boolean, is the gene to be DE?
#' @param model.zero.prob a boolean, should zero-density be modelled?
#' @param w weight
#' @param batch,group batch and grouping vectors 
#' @param prior.count the prior count for the cpm transofrm
#'
#' @return list of density estimates
geneParmEst <- function(cpm.data.i, batch, group, prior.count = prior.count,
                          de.ind, model.zero.prob, w){
  min.val = log(prior.count)
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
  n = length(Y)
  n.dist = length(unique(Y))
  if(n>=5 && n.dist>3){ # fit density only if there are at least 3 distinct values (e.g cpm)
    #Fit normal carrier density
    mu.hat = mean(Y)
    sig.hat = sd(Y)
    #Bin counts
    countY <- obtCount(Y = Y, w = w)
    if(sum(Y > min.val) > min.count.nonnull){
      #Fit exponential density
      fitLLmodel(yy = countY, mu.hat = mu.hat, sig.hat = sig.hat, n = n)
    } else{
      g0 = diff(pnorm(countY$breaks, mean = mu.hat, sd = sig.hat))
      g0[g0==0] = .Machine$double.eps
      c(countY, list(g0 = g0, n = n))
    }
  }else{
    c(0, list(g0 = .Machine$double.eps, n = n))
  }
  
}