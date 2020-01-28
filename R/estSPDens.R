# A function to estimate the semi-parametric densities (gene level)
estSPDens <- function(est.list.i, par.sample, DE.ind.ii, group, g0){
  if(DE.ind.ii==0){
    g1 <- lapply(seq_len(nrow(par.sample)), function(bb){ 
      b.data <- est.list.i$batch.est[[bb]]
      gg0 <- g0[[bb]]*sum(b.data$yy$Ny) +1
      
      s <- b.data$yy$S
      s.mat <- matrix(NA, ncol = length(coef(b.data$llm)), nrow=length(s))
      for(i in seq_len(ncol(s.mat))){
        s.mat[, i] <- s^(i-1)
      }
      gg1 <- exp(s.mat %*% par.sample[bb, seq_len(ncol(s.mat))])*gg0
      gg1 <- data.frame(gy=gg1, s=s, lls=b.data$yy$lls, uls=b.data$yy$uls)
      gg1$gy[is.infinite(gg1$gy)] <- max(gg1$gy[!is.infinite(gg1$gy)]) 
      gg1$Gy <- cumsum(gg1$gy)/sum(gg1$gy)
      gg1 
    })  
  }
  else{
    g1 <- lapply(sort(unique(group)), function(g){
      par.sample.g <- par.sample[[g]] #par.sample[[paste0("grp_", g)]]
      g0.g <- g0[[g]] #g0[[paste0("grp_", g)]]
      
      lapply(seq_len(nrow(par.sample.g)), function(bb){   
        b.data <- est.list.i$batch.est[[bb]][[g]]
        gg0 <- g0.g[[bb]]*sum(b.data$yy$Ny) + 1
        
        s <- b.data$yy$S
        s.mat <- matrix(NA, ncol = length(coef(b.data$llm)), nrow=length(s))
        for(k in seq_len(ncol(s.mat))){
          s.mat[, k] <- s^(k-1)
        }
        gg1 <- exp(s.mat %*% as.matrix(par.sample.g[bb, seq_len(ncol(s.mat))]))*gg0
        gg1 <- data.frame(gy=gg1, s=s, lls=b.data$yy$lls, uls=b.data$yy$uls)
        gg1$gy[is.infinite(gg1$gy)] <- max(gg1$gy[!is.infinite(gg1$gy)]) 
        gg1$Gy <- cumsum(gg1$gy)/sum(gg1$gy)
        gg1
      }) 
    }) 
  }
  return(g1)
}