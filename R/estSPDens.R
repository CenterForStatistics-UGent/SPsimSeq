# A function to estimate the semi-parametric densities (gene level)
estSPDens <- function(est.list.i, par.sample, DE.ind.ii, group, g0){
  if(DE.ind.ii==0){
    g1 <- lapply(seq_len(nrow(par.sample)), function(bb){ 
      b.data <- est.list.i[[bb]]
      gg0 <- g0[[bb]]*sum(b.data$Ny) +1
      s <- b.data$S
      s.mat <- buildXmat(s, nc = length(b.data$betas))
      gg1 <- exp(s.mat %*% par.sample[bb, seq_len(ncol(s.mat))])*gg0
      gg1 <- data.frame(gy=gg1, s=s, lls=b.data$lls, uls=b.data$uls)
      gg1$gy[is.infinite(gg1$gy)] <- max(gg1$gy[!is.infinite(gg1$gy)]) 
      gg1$Gy <- cumsum(gg1$gy)/sum(gg1$gy)
      gg1 
    })  
  }
  else{
    g1 <- lapply(sort(unique(group)), function(g){
      par.sample.g <- par.sample[[g]]
      g0.g <- g0[[g]]
      lapply(seq_len(nrow(par.sample.g)), function(bb){   
        b.data <- est.list.i[[bb]][[g]]
        gg0 <- g0.g[[bb]]*sum(b.data$Ny) + 1
        s <- b.data$S
        s.mat <- buildXmat(s, nc = length(b.data$betas))
        gg1 <- exp(s.mat %*% as.matrix(par.sample.g[bb, seq_len(ncol(s.mat))]))*gg0
        gg1 <- data.frame(gy=gg1, s=s, lls=b.data$lls, uls=b.data$uls)
        gg1$gy[is.infinite(gg1$gy)] <- max(gg1$gy[!is.infinite(gg1$gy)]) 
        gg1$Gy <- cumsum(gg1$gy)/sum(gg1$gy)
        gg1
      }) 
    }) 
  }
  return(g1)
}