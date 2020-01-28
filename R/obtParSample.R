# A function to estimate the parameters of the log-linear model (gene level)
obtParSample <- function(est.list.i, DE.ind.ii, n.batch, group){
  if(DE.ind.ii==0){
    if(length(n.batch)>1){   
      par.sample <- as.matrix(do.call("rbind.fill2", lapply(est.list.i$batch.est, function(bt){
        v.mat     <- bt$parm.list$v
        betas.vec <- bt$parm.list$betas 
        data.frame(t(as.matrix(c(mvrnorm(n = 1, mu= betas.vec, Sigma = v.mat), 
                                 mu.hat=bt$parm.list$mu.hat, sig.hat=bt$parm.list$sig.hat))))
      }))) 
    }
    else if(length(n.batch)==1){
      par.sample <- t(as.matrix(est.list.i$Mu.batch))
    } 
  }
  else{ 
    if(length(n.batch)>1){ 
      par.sample <- lapply(sort(unique(group)), function(g){
        
        par.sample.g <- as.matrix(do.call("rbind.fill2", 
                                          lapply(est.list.i$batch.est[[g]], function(bt){ 
                                            data.frame(t(as.matrix(c(mvrnorm(n = 1, 
                                                                             mu= bt$parm.list$betas, 
                                                                             Sigma = bt$parm.list$v), 
                                                                     mu.hat=bt$parm.list$mu.hat,
                                                                     sig.hat=bt$parm.list$sig.hat)))) 
                                          })))
        
        par.sample.g 
      }) 
    }
    else if(length(n.batch)==1){
      par.sample <- lapply(sort(unique(group)), function(g){
        par.sample.g <-t(as.matrix(est.list.i$Mu.batch[[g]]))
        
        par.sample.g 
      }) 
    } 
  }
  return(par.sample)
}