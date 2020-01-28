# A function to estimate the normal carrier density (gene level)
estCarrierDens <- function(est.list.i, par.sample, DE.ind.ii, group){
  if(DE.ind.ii==0){
    g0 <- lapply(seq_len(nrow(par.sample)), function(bb){
      b.data <- est.list.i$batch.est[[bb]] 
      gg0 <- (pnorm(b.data$yy$uls, as.matrix(par.sample)[bb, "mu.hat"],
                    as.matrix(par.sample)[bb, "sig.hat"]) -
                pnorm(b.data$yy$lls,as.matrix(par.sample)[bb, "mu.hat"], 
                      as.matrix(par.sample)[bb, "sig.hat"]))
      
      gg0[is.nan(gg0)] <- 0
      gg0
    })
  }
  else{
    g0 <- lapply(sort(unique(group)), function(g){
      par.sample.g <- par.sample[[g]] #par.sample[[paste0("grp_", g)]]
      lapply(seq_len(nrow(par.sample.g)), function(bb){  
        b.data <- est.list.i$batch.est[[bb]][[g]]
        gg0 <- (pnorm(b.data$yy$uls, par.sample.g[bb, "mu.hat"][[1]], 
                      par.sample.g[bb, "sig.hat"][[1]]) -
                  pnorm(b.data$yy$lls, par.sample.g[bb, "mu.hat"][[1]], 
                        par.sample.g[bb, "sig.hat"][[1]]))
        gg0[is.nan(gg0)] <- 0
        gg0 
      })
    }) 
  }
  return(g0)
}