# A function to estimate the parameters of the log-linear model (gene level)
obtParSample <- function(est.list.i, DE.ind.ii, n.batch, group){
  if(DE.ind.ii){
    if(length(n.batch)>1){   
      par.sample <- as.matrix(do.call("rbind", lapply(est.list.i, function(bt){
        v.mat     <- bt$v
        betas.vec <- bt$betas 
        data.frame(t(as.matrix(c(mvrnorm(n = 1, mu= betas.vec, Sigma = v.mat), 
                                 mu.hat=bt$mu.hat, sig.hat=bt$sig.hat))))
      }))) 
    }
    else if(length(n.batch)==1){
      par.sample <- t(as.matrix(c(est.list.i[[1]]$betas, mu.hat=est.list.i[[1]]$mu.hat, 
                                  sig.hat=est.list.i[[1]]$sig.hat)))
    } 
  }
  else{ 
    if(length(n.batch)>1){ 
      par.sample <- lapply(sort(unique(group)), function(g){
        est.bb <- lapply(est.list.i, function(bt){
          bt[[g]]
        })
        par.sample.g <- as.matrix(do.call("rbind", lapply(est.bb, function(btg){
          v.mat     <- btg$v
          betas.vec <- btg$betas 
          data.frame(t(as.matrix(c(mvrnorm(n = 1, mu= betas.vec, Sigma = v.mat), 
                                   mu.hat=btg$mu.hat, sig.hat=btg$sig.hat))))
        }))) 
        par.sample.g 
      })
      
      
    }
    else if(length(n.batch)==1){
      par.sample <- lapply(sort(unique(group)), function(g){
        t(as.matrix(c(est.list.i[[1]][[g]]$betas, mu.hat=est.list.i[[1]][[g]]$mu.hat, 
          sig.hat=est.list.i[[1]][[g]]$sig.hat)))
      }) 
    } 
  }
  return(par.sample)
}