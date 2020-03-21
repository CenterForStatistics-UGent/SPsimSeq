# A function to evaluate the normal carrier density (gene level)
evalCarrierDens <- function(est.list.i, par.sample, DE.ind.ii, group){
  g0 = if(!DE.ind.ii){
    g0 <- lapply(seq_len(nrow(par.sample)), function(bb){
      evalCarrierDensVec(est.list.i[[bb]], par.sample)
    })
  } else{
    lapply(sort(unique(group)), function(g){
      par.sample.g <- par.sample[[g]]
      lapply(seq_len(nrow(par.sample.g)), function(bb){  
        evalCarrierDensVec(b.data = est.list.i[[bb]][[g]], par.sample.g)
      })
    }) 
  }
  return(g0)
}
evalCarrierDensVec = function(b.data, par.sample){
  gg0 <- (pnorm(b.data$uls, par.sample[bb, "mu.hat"][[1]], 
                par.sample.g[bb, "sig.hat"][[1]]) -
            pnorm(b.data$lls, par.sample.g[bb, "mu.hat"][[1]], 
                  par.sample.g[bb, "sig.hat"][[1]]))
  gg0[is.nan(gg0)] <- 0
  gg0 
}
