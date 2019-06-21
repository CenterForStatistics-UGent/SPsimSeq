# Estimates the carrier density
SPestCarrier <- function(de.ind, par.sample, group, est.list.i, ...){
  if(de.ind==0){
    g0 <- lapply(1:nrow(par.sample), function(bb){
      b.data <- est.list.i$batch.est[[bb]] 
      gg0 <- (pnorm(b.data$yy$uls, as.matrix(par.sample)[bb, "mu.hat"],
                    as.matrix(par.sample)[bb, "sig.hat"]) -
                pnorm(b.data$yy$lls,as.matrix(par.sample)[bb, "mu.hat"],
                      as.matrix(par.sample)[bb, "sig.hat"]))
      # gg0 <- (pnorm(b.data$yy$Y+0.05, as.matrix(par.sample)[bb, "mu.hat"],
      #               as.matrix(par.sample)[bb, "sig.hat"]) -
      #           pnorm(b.data$yy$Y-0.05,as.matrix(par.sample)[bb, "mu.hat"], 
      #                 as.matrix(par.sample)[bb, "sig.hat"]))
      # 
      gg0[is.nan(gg0)] <- 0
      gg0
      # gg0.expd <- sapply(b.data$yy$Y, function(y) gg0[b.data$yy$lls<y & b.data$yy$uls>=y])
      # gg0.expd
    })
  }
  else{
    g0 <- lapply(sort(unique(group)), function(g){
      par.sample.g <- par.sample[[g]] #par.sample[[paste0("grp_", g)]]
      lapply(1:nrow(par.sample.g), function(bb){  
        b.data <- est.list.i$batch.est[[bb]][[g]] 
        gg0 <- (pnorm(b.data$yy$uls, par.sample.g[bb, "mu.hat"][[1]],
                      par.sample.g[bb, "sig.hat"][[1]]) -
                  pnorm(b.data$yy$lls, par.sample.g[bb, "mu.hat"][[1]],
                        par.sample.g[bb, "sig.hat"][[1]]))
        # gg0 <- (pnorm(b.data$yy$Y+0.05, par.sample.g[bb, "mu.hat"][[1]], 
        #               par.sample.g[bb, "sig.hat"][[1]]) -
        #           pnorm(b.data$yy$Y-0.05, par.sample.g[bb, "mu.hat"][[1]], 
        #                 par.sample.g[bb, "sig.hat"][[1]]))
        gg0[is.nan(gg0)] <- 0
        gg0
        # gg0.expd <- sapply(b.data$yy$Y, function(y) gg0[b.data$yy$lls<y & b.data$yy$uls>=y])
        # gg0.expd
      })
    })
    #names(g0) <- paste0("grp_", sort(unique(group)))
  }
  g0
}