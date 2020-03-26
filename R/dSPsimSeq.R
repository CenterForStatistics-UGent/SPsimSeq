#' Evaluate the densities in the estimated SPsimSeq object 
#' @param SPobj The SPsimSeq object, with details retained
#' @param newData the data points at which the density should be evaluated
#'@param force.fit.data a logical value. If TRUE (the default), then observations larger or smaller
#'than the midpoint of the higher or lower classes (respectively) will be forced to the midpoint of
#'the higher or lower classes (respectively). This argument is particularly important to avoid
#'extrapolation problem for observations that are out of the range of the data points (midpoints)
#'used for the estimation of the SPsimSeq log-linear regression parameters.
#'@return a list of estimated densities
#'@export
#'@examples
#' data("zhang.data.sub")
#' # filter genes with sufficient expression (important step to avoid bugs)
#' zhang.counts <- zhang.data.sub$counts[rowSums(zhang.data.sub$counts > 0)>=5, ]
#' MYCN.status  <- zhang.data.sub$MYCN.status
#' # simulate data
#' zhang.counts2 <- zhang.counts[sample(nrow(zhang.counts), 2000), ]
#' sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = zhang.counts2,
#'                           group = MYCN.status, n.genes = 2000, batch.config = 1,
#'                           group.config = c(0.5, 0.5), tot.samples = 20,
#'                           pDE = 0.1, lfc.thrld = 0.5, result.format = "list",
#'                           return)
#' outDens = dSPsimSeq(sim.data.bulk, newData = zhang.counts2[1,])
# select.genes <- sample(rownames(cpm.dat), 20)
# #win.graph()
# par(mfrow=c(4, 5))
# for(i in select.genes){
#   if(length(cpm.dat.f[[i]])==1){
#     plot(cpm.dat[i, ], cpm.dat.f[[i]][[1]], type = "p", main = i)
#     lines(density(cpm.dat[i, ]), col="gray")
#   }else{
#     plot(cpm.dat[i, ], cpm.dat.f[[i]][[1]], type = "p", ylim = c(min(min(cpm.dat.f[[i]][[1]]),
#                                                                      min(cpm.dat.f[[i]][[2]])),
#                                                                  max(max(cpm.dat.f[[i]][[1]]),
#                                                                      max(cpm.dat.f[[i]][[2]]))))
#     points(cpm.dat[i, ], cpm.dat.f[[i]][[2]], type = "p", col=2)
#   }
# }
dSPsimSeq <- function(SPobj, est.parms, force.fit.data=TRUE){
  if(is.list(est.parms$Mu.batch)){
    n.groups <- length(est.parms$Mu.batch) 
  }else{
    n.groups <- 1
  }
  n.groups <- length(est.parms[[1]]) 
  
  if(n.groups==1){ 
    if(force.fit.data){
      x[x<min(est.parms$batch.est[[1]]$S)] <- min(est.parms$batch.est[[1]]$S)
      x[x>max(est.parms$batch.est[[1]]$S)] <- max(est.parms$batch.est[[1]]$S)
    }
    g0 <- dnorm(x, est.parms$Mu.batch[["mu.hat"]], 
                 est.parms$Mu.batch[["sig.hat"]])
    beta.hat <- est.parms$Mu.batch[!(names(est.parms$Mu.batch) %in% 
                                                 c("mu.hat", "sig.hat"))]
    x.mat <- buildXmat(x, nc = length(beta.hat))
    g1 <- exp(x.mat%*%as.matrix(beta.hat)) 
    f.hat <- list((g0)*g1)
  }else if(n.groups>1){
    f.hat <- lapply(seq_len(n.groups), function(ii){
      if(force.fit.data){
        x[x<min(est.parms[[1]][[ii]]$S)] <- 
          min(est.parms[[1]][[ii]]$S)
        x[x>max(est.parms[[1]][[ii]]$S)] <- 
          max(est.parms[[1]][[ii]]$S)
      }
      g0 <- dnorm(x, est.parms[[1]][[ii]]$mu.hat, 
                  est.parms[[1]][[ii]]$sig.hat) 
      beta.hat <- est.parms[[1]][[ii]]$betas
      x.mat <- buildXmat(x, nc = length(beta.hat))
      g1 <- exp(x.mat%*%as.matrix(beta.hat))
      g0*g1 
    })
    names(f.hat) <- names(est.parms$SPsim.dens.parms)
  }
  f.hat
}