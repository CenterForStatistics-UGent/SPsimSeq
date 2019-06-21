SPestDens <- function(de.ind, group, batch, par.sample, est.list.i, g0, ...){
  
  if(de.ind==0){
    g1 <- lapply(1:nrow(par.sample), function(bb){ 
      b.data <- est.list.i$batch.est[[bb]]
      gg0 <- g0[[bb]]*sum(b.data$yy$Ny) +1
      SPestDensBB(b.data = b.data, gg0 = gg0, par.sample = par.sample, bb = bb)
      
      # s <- b.data$yy$Y
      # s.mat <- matrix(NA, ncol = length(coef(b.data$llm)), nrow=length(s))
      # for(ii in 1:ncol(s.mat)){
      #   s.mat[, ii] <- s^(ii-1)
      # }
      # gg1 <- exp(s.mat %*% par.sample[bb, 1:ncol(s.mat)])*gg0
      # #gg1 <- data.frame(gy=gg1, s=s, lls=b.data$yy$lls, uls=b.data$yy$uls)
      # lls.expd <- sapply(s, function(ss) b.data$yy$lls[b.data$yy$lls<ss & b.data$yy$uls>=ss])
      # uls.expd <- sapply(s, function(ss) b.data$yy$uls[b.data$yy$lls<ss & b.data$yy$uls>=ss])
      # S.expd   <- sapply(s, function(ss) b.data$yy$S[b.data$yy$lls<ss & b.data$yy$uls>=ss])
      # Ny.expd  <- sapply(s, function(ss) b.data$yy$Ny[b.data$yy$lls<ss & b.data$yy$uls>=ss])
      # gg1 <- data.frame(s=s, S.expd=S.expd, lls.expd=lls.expd, 
      #                   uls.expd=uls.expd, Ny.expd=Ny.expd, gg0=gg0, gy=gg1)
      # gg1 <- gg1[!duplicated(gg1),]
      # gg1$gy[is.infinite(gg1$gy)] <- max(gg1$gy[!is.infinite(gg1$gy)]) 
      # gg1 <- gg1[order(gg1$s),]
      # gg1$Gy <- cumsum(gg1$gy/gg1$Ny.expd)/sum(gg1$gy/gg1$Ny.expd) 
      # gg1
      
      
    })  
    # plot(g1[[1]]$s, g1[[1]]$gy, type="b", xlim=range(do.call('c', lapply(g1, function(x) x[, "s"]))))
    # for(k in 2:length(g1)){
    #   lines(g1[[k]]$s, g1[[k]]$gy, type="b", col=k)
    # }
  }
  else{
    g1 <- lapply(sort(unique(group)), function(g){
      par.sample.g <- par.sample[[g]] #par.sample[[paste0("grp_", g)]]
      g0.g <- g0[[g]] #g0[[paste0("grp_", g)]]
      
      lapply(1:nrow(par.sample.g), function(bb){   
        b.data <- est.list.i$batch.est[[bb]][[g]]
        gg0 <- g0.g[[bb]]*sum(b.data$yy$Ny) + 1
        
        SPestDensBB(b.data = b.data, gg0 = gg0, par.sample = par.sample.g, bb = bb)
        # s <- b.data$yy$Y 
        # s.mat <- matrix(NA, ncol = length(coef(b.data$llm)), nrow=length(s))
        # for(k in 1:ncol(s.mat)){
        #   s.mat[, k] <- s^(k-1)
        # }
        # gg1 <- exp(s.mat %*% as.matrix(par.sample.g[bb, 1:ncol(s.mat)]))*gg0 
        # #gg1 <- data.frame(gy=gg1, s=s, lls=b.data$yy$lls, uls=b.data$yy$uls)
        # lls.expd <- sapply(s, function(ss) b.data$yy$lls[b.data$yy$lls<ss & b.data$yy$uls>=ss])
        # uls.expd <- sapply(s, function(ss) b.data$yy$uls[b.data$yy$lls<ss & b.data$yy$uls>=ss])
        # S.expd   <- sapply(s, function(ss) b.data$yy$S[b.data$yy$lls<ss & b.data$yy$uls>=ss])
        # Ny.expd  <- sapply(s, function(ss) b.data$yy$Ny[b.data$yy$lls<ss & b.data$yy$uls>=ss])
        # 
        # gg1 <- data.frame(s=s, S.expd=S.expd, lls.expd=lls.expd, 
        #                   uls.expd=uls.expd, Ny.expd=Ny.expd, gg0=gg0, gy=gg1)
        # gg1 <- gg1[!duplicated(gg1),]
        # gg1$gy[is.infinite(gg1$gy)] <- max(gg1$gy[!is.infinite(gg1$gy)]) 
        # gg1 <- gg1[order(gg1$s),]
        # gg1$Gy <- cumsum(gg1$gy/gg1$Ny.expd)/sum(gg1$gy/gg1$Ny.expd) 
        # gg1
      }) 
    })
    #names(g1) <- paste0("grp_", sort(unique(group)))
  }
  g1
}