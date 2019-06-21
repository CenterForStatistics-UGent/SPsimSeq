# Estimates the density in each batch for a gene i
SPestDensBB <- function(b.data, gg0, par.sample, bb){ 
  s <- b.data$yy$S
  
  s.mat <- matrix(NA, ncol = length(coef(b.data$llm)), nrow=length(s))
  for(ii in 1:ncol(s.mat)){
    s.mat[, ii] <- s^(ii-1)
  }
  gg1 <- exp(s.mat %*% par.sample[bb, 1:ncol(s.mat)])*gg0
  gg1 <- data.frame(gy=gg1, s=s, lls=b.data$yy$lls, uls=b.data$yy$uls, Ny=b.data$yy$Ny)
  # lls.expd <- sapply(s, function(ss) b.data$yy$lls[b.data$yy$lls<ss & b.data$yy$uls>=ss])
  # uls.expd <- sapply(s, function(ss) b.data$yy$uls[b.data$yy$lls<ss & b.data$yy$uls>=ss])
  # S.expd   <- sapply(s, function(ss) b.data$yy$S[b.data$yy$lls<ss & b.data$yy$uls>=ss])
  # Ny.expd  <- sapply(s, function(ss) b.data$yy$Ny[b.data$yy$lls<ss & b.data$yy$uls>=ss])
  # gg1 <- data.frame(s=s, S.expd=S.expd, lls.expd=lls.expd,
  #                   uls.expd=uls.expd, Ny.expd=Ny.expd, gg0=gg0, gy=gg1)
  #  
  gg1 <- gg1[!duplicated(gg1),]
  gg1$gy[is.infinite(gg1$gy)] <- max(gg1$gy[!is.infinite(gg1$gy)]) 
  gg1 <- gg1[order(gg1$s),]
  gg1$Gy <- cumsum(gg1$gy/gg1$Ny)/sum(gg1$gy/gg1$Ny) 
  gg1
}