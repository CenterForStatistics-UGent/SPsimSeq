
## This function is suggested only for a particular problem, and not to be used in general
#Density estimation based on SPsimSeq parameters 
#
#@param x a numerical vector of observations for which we want to estimate the density
#@param est.parms a list object returned from SPsimSeq function with name 'SPsim.est.densities'
#for a particular gene/feature (see the example for details)
#@param force.fit.data a logical value. If TRUE (the default), then observations larger or smaller
#than the midpoint of the higher or lower classes (respectively) will be forced to the midpoint of
#the higher or lower classes (respectively). This argument is particularly important to avoid
#extrapolation problem for observations that are out of the range of the data points (midpoints)
#used for the estimation of the SPsimSeq log-linear regression parameters.
#
#@return a list of estimated densities
#
#
#@examples (below)
###
dSPsimSeq <- function(x, est.parms, force.fit.data=TRUE){
  n.group <- length(est.parms$SPsim.dens.hat)
  if(n.group==1){ 
    if(force.fit.data){
      x[x<min(est.parms$SPsim.dens.hat[[1]]$s)] <- min(est.parms$SPsim.dens.hat[[1]]$s)
      x[x>max(est.parms$SPsim.dens.hat[[1]]$s)] <- max(est.parms$SPsim.dens.hat[[1]]$s)
    }
    g0 <- dnorm(x, est.parms$SPsim.dens.parms[,"mu.hat"], 
                 est.parms$SPsim.dens.parms[,"sig.hat"])
    beta.hat <- est.parms$SPsim.dens.parms[, !(colnames(est.parms$SPsim.dens.parms) %in% 
                                                 c("mu.hat", "sig.hat"))]
    x.mat <- buildXmat(x, nc = length(beta.hat))
    g1 <- exp(x.mat%*%as.matrix(beta.hat)) 
    f.hat <- list((g0)*g1)
  }else if(n.group>1){
    f.hat <- lapply(seq_len(n.group), function(ii){
      if(force.fit.data){
        x[x<min(est.parms$SPsim.dens.hat[[ii]][[1]]$s)] <- 
          min(est.parms$SPsim.dens.hat[[ii]][[1]]$s)
        x[x>max(est.parms$SPsim.dens.hat[[ii]][[1]]$s)] <- 
          max(est.parms$SPsim.dens.hat[[ii]][[1]]$s)
      }
      g0 <- dnorm(x, est.parms$SPsim.dens.parms[[ii]][,"mu.hat"], 
                  est.parms$SPsim.dens.parms[[ii]][,"sig.hat"]) 
      beta.hat <- est.parms$SPsim.dens.parms[[ii]][, 
                      !(colnames(est.parms$SPsim.dens.parms[[ii]]) %in% c("mu.hat", "sig.hat"))]
      x.mat <- buildXmat(x, nc = length(beta.hat))
      g1 <- exp(x.mat%*%as.matrix(beta.hat))
      g0*g1 
    })
    names(f.hat) <- names(est.parms$SPsim.dens.parms)
  }
  f.hat
}
# Example
# load(".../problemData.RData")
# mcr_data <- mat
# mcr_data <- mcr_data[rownames(mcr_data) != "", ]
# 
# sim.data <- SPsimSeq(n.sim = 1, s.data = mcr_data, batch = NULL, genewiseCor = FALSE,
#                           group = fac, n.genes = 4900, batch.config = 1,
#                           group.config = prop.table(table(fac)), tot.samples = 30,
#                           pDE = 0.2, lfc.thrld = 0.5, t.thrld = 2.0, llStat.thrld = 0,
#                           w=0.6, log.CPM.transform = FALSE,
#                           result.format = "list",  seed = 2581988)
# 
# sim.data1 <- sim.data[[1]]
# 
# cpm.dat <- mcr_data
# cpm.dat <- cpm.dat[as.character(unique(sim.data1$rowData$source.ID)), ]
# 
# cpm.dat.f <- lapply(rownames(cpm.dat), function(rr){
#   #print(rr)
#   rr2 <- rownames(sim.data1$rowData)[sim.data1$rowData$source.ID==rr]
#   rr2 <- rr2[1]
#   dSPsimSeq(x = cpm.dat[rr, ], est.parms = sim.data1$SPsim.est.densities[[rr2]])
# })
# names(cpm.dat.f) <- rownames(cpm.dat)
# 
# select.genes <- sample(rownames(cpm.dat), 20)
# #win.graph()
# par(mfrow=c(4, 5))
# for(i in select.genes){
#   if(length(cpm.dat.f[[i]])==1){
#     plot(cpm.dat[i, ], cpm.dat.f[[i]][[1]], type = "p", main = i)
#     #lines(density(cpm.dat[i, ]), col="gray")
#   }else{
#     plot(cpm.dat[i, ], cpm.dat.f[[i]][[1]], type = "p", ylim = c(min(min(cpm.dat.f[[i]][[1]]),
#                                                                      min(cpm.dat.f[[i]][[2]])),
#                                                                  max(max(cpm.dat.f[[i]][[1]]),
#                                                                      max(cpm.dat.f[[i]][[2]]))))
#     points(cpm.dat[i, ], cpm.dat.f[[i]][[2]], type = "p", col=2)
#   }
# }
