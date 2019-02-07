#' A function to simulate bulk or single cell RNA sequencing data
#' 
#' @param n.sim number of simulated data to be created
#' @param s.data a source data (a SingleCellExperiment object)
#' @param batch a vector containg btach indicator for each sample/cell
#' @param group a vector containg group indicator for each sample/cell
#' @param cand.genes a list object contating candidate genes 
#' @param n.genes a numeric value for the total number of genes to be simulated
#' @param pDE a numeric value between 0 and 1 indicating the fraction of DE genes 
#' in a single simulated data
#' @param batch.config a numerical vector for the marginal fraction of samples in each batch. 
#' The number of batches to be simulated is equal to the size of the vector.
#' All values must sum to 1.
#' @param group.config a numerical vector for the marginal fraction of samples in each group.
#' The number of groups to be simulated is equal to the size of the vector. All values must sum to 1.
#' @param model.zero.prob a logical value whether to model the zero probablity separately 
#' (suitable for single cell data)
#' @param tot.samples total number of samples to be simulated.
#' @param null.group a character value that indicating a name of
#' @param fc.type a character indicating how DE genes should be calculated 
#' ('mean.diff'= difference in the mean log CPM, 'rank.sum'= a U-statistic for the log CPM)
#' @param lfc.thrld a numeric value for the minimum fold change for DE genes (if fc.type='mean.diff')
#' @param t.thrld a numeric value for the minimum t statistic for DE genes (if fc.type='mean.diff')
#' @param U.thrld a numeric value for the minimum U-statistic for DE genes (if fc.type='rank.sum')
#' @param llStat.thrld a numeric value for the minimum squared test statistics from a log-linear 
#' model containing X as a covariate to select DE genes
#' @param carrier.dist a character indicating the type of carrier density (carrier.dist="normal" 
#' or carrier.dist="kernel")
#' @param w a numeric value between 0 and 1 or NULL refering the number of classes to be created 
#' for the outcome data (if NULL the algorithm to calculate breakes in graphics::hist() function 
#' will be used)
#' @param max.frac.zero a numeric value between 0 and 1 indicating the maximum fraction of zero 
#' counts that a DE gene should have
#' @param max.frac.zeror.diff a numeric value between 0 and 1 indicating the maximum  absolute 
#' difference in the fraction of zero counts between the groups for DE genes
#' @param verbose a logical value, if TRUE it displays a message about the satatus of the simulation
#' @param n.mean.class an integer for the number of groups to be created for the mean log CPM of genes
#' @param subset.data a logical value to subset columns(samples) if the size of the data is too big
#' for space and computation time saving
#' @param n.samples an integer indicating the number of samples/cells to sample (if subset.data=TRUE)
#' @param  ... further arguments passed to or from other methods.
#' 
#' @return a list of SingleCellExperiment object  each contatining simulated counts (not normalized), 
#' cell level information in colData, and gene level information in rowData
#' 
#' @examples 
#' ##
#' @export 
#' 
SPsimSeqData <- function(n.sim=1,   s.data, batch=NULL, group=NULL, cand.genes=NULL, n.genes=1000,  
                       batch.config=c(0.33, 0.33, 0.34),  group.config=c(0.5, 0.5), 
                       tot.samples=150, pDE=0.2, model.zero.prob=FALSE,
                       null.group=which.max(table(group)), 
                       fc.type="mean.diff",  lfc.thrld=0.5, U.thrld=0.7, llStat.thrld=5, 
                       t.thrld=2.5,  carrier.dist="normal", w=0.5, max.frac.zero=0.7, 
                       max.frac.zeror.diff=0.1, n.mean.class=30, 
                       subset.data=FALSE, n.samples=400, verbose=TRUE, ...)
{
  
  # experiment configurartion
  if(verbose) {message("Configuring design ...")}
  exprmt.design <- expriment.config(batch.config = batch.config, group.config = group.config,
                                    tot.samples = tot.samples)
  n.batch <- exprmt.design$n.batch
  n.group <- exprmt.design$n.group
  config.mat <- exprmt.design$exprmt.config
  
  # prepare source data
  if(verbose) {message("Preparing source data ...")}
  prepare.S.Data <- prepareSourceData(s.data, batch = batch, group = group, cand.genes = cand.genes,
                                      exprmt.design=exprmt.design, fc.type=fc.type, lfc.thrld=lfc.thrld,
                                      U.thrld=U.thrld, llStat.thrld=llStat.thrld, t.thrld=t.thrld,
                                      carrier.dist=carrier.dist, w=w, max.frac.zero=max.frac.zero)
  LL <- prepare.S.Data$LL
  cpm.data <- prepare.S.Data$cpm.data
  sub.batchs <- prepare.S.Data$sub.batchs
  
  
  # fit logistic regression for the probability of zeros
  if(model.zero.prob){
    if(verbose) {message("Fitting zero probability model ...")}
    if(class(s.data)=="SingleCellExperiment"){
      fracZero.logit.list <- lapply(unique(sub.batchs), function(b){
        zeroProbModel(counts(s.data)[, batch==b], subset.data=subset.data, n.samples=n.samples)
      }) 
    }
    else if(class(s.data) %in% c("data.frame", "matrix")){
      fracZero.logit.list <- lapply(unique(sub.batchs), function(b){
        zeroProbModel(s.data[, batch==b], subset.data=subset.data, n.samples=n.samples)
      })
    }
  }
  
  
  # candidate genes
  if(verbose) {message("Selecting genes ...")}
  null.genes0     <- prepare.S.Data$cand.genes$null.genes
  nonnull.genes0  <- prepare.S.Data$cand.genes$nonnull.genes
  
  if(pDE>0 & !is.null(group) & length(group.config)>1){
    if((1-pDE)*n.genes <= length(null.genes0)){
      null.genes     <- sample(null.genes0, (1-pDE)*(n.genes), replace = FALSE)
    }
    else{
      null.genes     <- sample(null.genes0, (1-pDE)*(n.genes), replace = TRUE)
    }
    
    if(pDE*n.genes <= length(nonnull.genes0)){
      nonnull.genes  <- sample(nonnull.genes0, pDE*n.genes, replace = FALSE)
    }
    else{
      nonnull.genes  <- sample(nonnull.genes0, pDE*n.genes, replace = TRUE)
    }
    
    sel.genes <- c(null.genes, nonnull.genes)
    DE.ind <- ifelse(sel.genes %in% null.genes, 0, 1)
    names(DE.ind) <- sel.genes
  }
  else{
    if(n.genes <= length(null.genes0)){
      sel.genes     <- sample(null.genes0, n.genes, replace = FALSE)
    }
    else{
      sel.genes     <- sample(null.genes0, n.genes, replace = TRUE)
    } 
    DE.ind <- rep(0, length(sel.genes))
    names(DE.ind) <- sel.genes
  } 
  
  # simulation step
  if(verbose) {message("Simulating data ...")}
  sim.data.list <- lapply(1:n.sim, function(h){
    if(verbose) {message(" ...", h, "of", n.sim)}
    # estimate batch specific parameters
    est.list <- lapply(sel.genes, function(i){
      #print(i)
      batch.est <- lapply(unique(batch), function(b){
        if(DE.ind[i]==0){
          Y0 <- cpm.data[i, (batch==b & group==null.group)] 
          if(model.zero.prob) Y <- Y0[Y0>0]
          else Y <- Y0
          
          countY   <- obtCount(Y)
          # plot(countY$S, countY$Ny) 
          parm.est <- fitLLmodel(countY)
          if(!is.null(parm.est$parm.list$betas)){
            parm.est
          } 
          else{NULL}
        }
        else{
          Y0 <- split(cpm.data[i, batch==b], group[batch==b])
          Y <- lapply(Y0, function(y0){
            if(model.zero.prob) {y0[y0>0]}
            else {y0}
          })
          
          countY   <- lapply(Y,  obtCount)
          # plot(countY$`1`$S, countY$`1`$Ny, type="l", xlim=range(cpm.data[i, batch==b])) 
          # lines(countY$`2`$S, countY$`2`$Ny, col=2) 
          parm.est <- lapply(countY, fitLLmodel)
          if(all(sapply(parm.est, function(x) !is.null(x$parm.list$betas)))){
            parm.est
          } 
          else{NULL}
        }
      }) 
      
      if(DE.ind[i]==0){
        batch.parms <- t(sapply(batch.est, function(b){ 
          if(!is.null(b)){
            c(b$parm.list$betas, mu.hat=b$parm.list$mu.hat, 
              sig.hat=b$parm.list$sig.hat)
          }
          else{NULL}
        }))
        if(!is.null(batch.est)){
          Mu.batch <- colMeans(batch.parms)
          V.batch  <- var(batch.parms) 
        }
        else{
          Mu.batch <- NULL
          V.batch  <- NULL
        }
      }
      else{
        batch.parms <- lapply(unique(group), function(g){
          b <- lapply(batch.est, function(x) x[[g]])
          b <- t(sapply(b, function(bb){
            if(!is.null(bb)){
              c(bb$parm.list$betas, mu.hat=bb$parm.list$mu.hat, 
                sig.hat=bb$parm.list$sig.hat)
            }
            else{NULL}
          })) })
        if(!any(sapply(batch.parms, is.null))){
          Mu.batch <- lapply(batch.parms, colMeans)
          V.batch  <- lapply(batch.parms, var)
        }
        else{
          Mu.batch <- NULL
          V.batch  <- NULL
        }
      }
      
      if(!is.null(Mu.batch)){
        list(Mu.batch=Mu.batch, V.batch=V.batch,
             batch.est=batch.est)
      }
      else{
        NULL
      }
      
    })
    
    sim.list <- lapply(1:length(sel.genes), function(i){
      if(!is.null(est.list[[i]])){
        #print(i)
        # sample params from MVN 
        if(DE.ind[i]==0){
          if(length(n.batch)>1){
            par.sample <- mvrnorm(n = length(sub.batchs), mu= est.list[[i]]$Mu.batch,
                                  Sigma = est.list[[i]]$V.batch) 
          }
          else if(length(n.batch)==1){
            par.sample <- t(as.matrix(est.list[[i]]$Mu.batch))
          } 
        }
        else{ 
          if(length(n.batch)>1){
            par.sample <- lapply(unique(group), function(g){
              mvrnorm(n = length(sub.batchs), mu= est.list[[i]]$Mu.batch[[g]],
                      Sigma = est.list[[i]]$V.batch[[g]]) 
            })
          }
          else if(length(n.batch)==1){
            par.sample <- lapply(unique(group), function(g){
              t(as.matrix(est.list[[i]]$Mu.batch[[g]]))
            })
          } 
        }
        
        
        # estimate carrier density (g0)
        if(DE.ind[i]==0){
          g0 <- lapply(1:nrow(par.sample), function(bb){
            b.data <- est.list[[i]]$batch.est[[bb]]
            gg0 <- (pnorm(b.data$yy$uls, par.sample[bb, "mu.hat"][[1]], par.sample[bb, "sig.hat"][[1]]) -
                      pnorm(b.data$yy$lls, par.sample[bb, "mu.hat"][[1]], par.sample[bb, "sig.hat"][[1]]))
            gg0[is.nan(gg0)] <- 0
            gg0
          })
        }
        else{
          g0 <- lapply(unique(group), function(g){
            par.sample.g <- par.sample[[g]]
            lapply(1:nrow(par.sample.g), function(bb){ 
              b.data <- est.list[[i]]$batch.est[[bb]][[g]]
              gg0 <- (pnorm(b.data$yy$uls, par.sample.g[bb, "mu.hat"][[1]], par.sample.g[bb, "sig.hat"][[1]]) -
                        pnorm(b.data$yy$lls, par.sample.g[bb, "mu.hat"][[1]], par.sample.g[bb, "sig.hat"][[1]]))
              gg0[is.nan(gg0)] <- 0
              gg0
            })
          })
        }
        
        
        
        # estimate the density g1(y)
        if(DE.ind[i]==0){
          g1 <- lapply(1:nrow(par.sample), function(bb){
            b.data <- est.list[[i]]$batch.est[[bb]]
            gg0 <- g0[[bb]]*sum(b.data$yy$Ny) +1
            
            s <- b.data$yy$S
            s.mat <- matrix(NA, ncol = length(coef(b.data$llm)), nrow=length(s))
            for(i in 1:ncol(s.mat)){
              s.mat[, i] <- s^(i-1)
            }
            gg1 <- exp(s.mat %*% as.matrix(par.sample[bb, 1:ncol(s.mat)]))*gg0
            gg1 <- data.frame(gy=gg1, s=s, lls=b.data$yy$lls, uls=b.data$yy$uls)
            gg1$gy[is.infinite(gg1$gy)] <- max(gg1$gy[!is.infinite(gg1$gy)]) 
            gg1$Gy <- cumsum(gg1$gy)/sum(gg1$gy)
            gg1 
          })  
          # plot(g1[[1]]$s, g1[[1]]$gy, type="b", xlim=range(do.call('c', lapply(g1, function(x) x[, "s"]))))
          # for(k in 2:length(g1)){
          #   lines(g1[[k]]$s, g1[[k]]$gy, type="b", col=k)
          # }
        }
        else{
          g1 <- lapply(unique(group), function(g){
            par.sample.g <- par.sample[[g]]
            g0.g <- g0[[g]]
            
            lapply(1:nrow(par.sample.g), function(bb){
              b.data <- est.list[[i]]$batch.est[[bb]][[g]]
              gg0 <- g0.g[[bb]]*sum(b.data$yy$Ny) +1
              
              s <- b.data$yy$S
              s.mat <- matrix(NA, ncol = length(coef(b.data$llm)), nrow=length(s))
              for(k in 1:ncol(s.mat)){
                s.mat[, k] <- s^(k-1)
              }
              gg1 <- exp(s.mat %*% as.matrix(par.sample.g[bb, 1:ncol(s.mat)]))*gg0
              gg1 <- data.frame(gy=gg1, s=s, lls=b.data$yy$lls, uls=b.data$yy$uls)
              gg1$gy[is.infinite(gg1$gy)] <- max(gg1$gy[!is.infinite(gg1$gy)]) 
              gg1$Gy <- cumsum(gg1$gy)/sum(gg1$gy)
              gg1 
            }) 
          })
        }
        
        # simulate data 
        if(DE.ind[i]==0){
          Y.star <- lapply(1:nrow(par.sample), function(bb){
            u <- runif(n.batch[[bb]])
            gg1 <- g1[[bb]]
            y.star.b <- sapply(u, function(uu){
              yy  <- gg1$s[which.min(abs(gg1$Gy-uu))]
              eps <- gg1$s[2]-gg1$s[1] #=mean(diff(gg1$s))
              yy  <- runif(1, yy-eps/2, yy+eps/2)
              yy
            })
            LL.b <- as.numeric(do.call("c", LL[[bb]]))
            y.star2.b <- round(((exp(y.star.b)-1)*LL.b)/1e6)
            y.star2.b[y.star2.b<0] <- 0
            #y.star2.b
            
            if(model.zero.prob){
              lLL_b <- log(LL.b)
              pred.pz  <- predict(fracZero.logit.list[[bb]], type="response",
                                  newdata=data.frame(x1=mean(y.star.b), x2=lLL_b))
              drop.mlt <- sapply(pred.pz, function(p) rbinom(1, 1, p))
              y.star2.b <- y.star2.b*(1-drop.mlt)
            }
            as.numeric(y.star2.b)
          })
          Y.star <- do.call("c", Y.star)
          #if(any(is.na(Y.star))){print(i)}
        }
        else{
          Y.star <- lapply(unique(group), function(g){
            par.sample.g <- par.sample[[g]]
            g1.g <- g1[[g]]
            y.star.g <- lapply(1:nrow(par.sample.g), function(bb){
              u <- runif(config.mat[bb, g])#runif(n.batch[[bb]]/length(n.group))
              gg1 <- g1.g[[bb]]
              y.star.b <- sapply(u, function(uu){
                yy  <- gg1$s[which.min(abs(gg1$Gy-uu))]
                eps <- gg1$s[2]-gg1$s[1] #=mean(diff(gg1$s))
                yy  <- runif(1, yy-eps/2, yy+eps/2)
                yy
              })
              LL.b.g <- LL[[bb]][[g]]
              y.star2.b <- round(((exp(y.star.b)-1)*LL.b.g)/1e6)
              y.star2.b[y.star2.b<0] <- 0
              
              if(model.zero.prob){
                lLL.b.g<- log(LL.b.g)
                pred.pz  <- predict(fracZero.logit.list[[bb]], type="response",
                                    newdata=data.frame(x1=mean(y.star.b), x2=lLL.b.g))
                drop.mlt <- sapply(pred.pz, function(p) rbinom(1, 1, p))
                y.star2.b <- y.star2.b*(1-drop.mlt)
              }
              as.numeric(y.star2.b)
            })
            y.star.g <- do.call("c", y.star.g)
          })
          Y.star <- do.call("c", Y.star)
          #if(any(is.na(Y.star))){print(i)}
        }}
      else{
        #print(i)
        Y.star <- rep(0, tot.samples)
        Y.star
      }
    }) 
    sim.count <- do.call(rbind, sim.list)
    rownames(sim.count) <- paste0("Gene_", 1:nrow(sim.count))
    colnames(sim.count) <- paste0("Sample_", 1:ncol(sim.count))
    
    col.data <- data.frame(Batch=rep(rep(1:length(n.batch), times=length(n.group)), config.mat),
                           Group=rep(rep(1:length(n.group), each=length(n.batch)), config.mat), 
                           sim.Lib.Size=do.call("c", do.call("c", LL)),   LS = colSums(sim.count),
                           row.names = colnames(sim.count))
    row.data <- data.frame(DE.ind=DE.ind, row.names = rownames(sim.count))
    sim.data <- SingleCellExperiment(assays=list(counts=sim.count),
                                     colData=col.data, rowData=row.data)
    sim.data
  })
  sim.data.list
}