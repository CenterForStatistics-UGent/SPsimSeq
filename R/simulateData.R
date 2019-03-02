#' A function to simulate bulk or single cell RNA sequencing data
#' 
#' @description This function simulates RNA sequencing data given a real RNA-seq data using
#' semi-parametric density estimation.
#' 
#' @param n.sim a numerical value for the number of simulated data to be  generated
#' @param s.data a source data (a SingleCellExperiment class object or a matrix/data.frame of counts with genes in 
#' rows and samples in columns)
#' @param batch a vector containg btach indicator for each sample/cell
#' @param group a vector containg group indicator for each sample/cell 
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
#' @param tot.samples a numerical value for total number of samples to be simulated. 
#' @param result.format a character value for the type of format for the output. Choice can  be
#' 'SCE' for SingleCellExperiment class or "list" for a list object that contains the simulated count,
#' column information abd row information.
#' @param verbose a logical value, if TRUE it displays a message about the satatus of the simulation
#' @param  ... further arguments passed to or from other methods.
#' 
#' @return a list of SingleCellExperiment object each contatining simulated counts (not normalized), 
#' cell level information in colData, and gene level information in rowData.
#' 
#' 
#' @details This function constructs the actual density of a given bulk or single cell RNA-seq data 
#' (passed using \emph{s.data} argument) using a specially designed exponetial family for density
#' estimation. Afterwards, a new data set is simulated from the estimated density. The objective is 
#' to maximaly retainthe characterisics of the source data in the simulated data without distributional 
#' assumption to the gene expression. This is particulaly useful for cases with complex/unknown 
#' distributions of the data, such as for single cell RNA-seq.
#' 
#' For single cell RNA-seq data, the high abundance of zero counts is a special characteristic as
#' a result of either biological or technical factors. Therefore, this package has a special feature to
#' model the probability of zero counts given the mean expression and library size  to add excess zeros
#' representing the technical effectss. This can be achieved by using \emph{model.zero.prob=TRUE}. Note that,
#' for extremly large size data, this feature may not work and a rrandom sample of cells (up to 400) can
#' be used to fit the model. To enable sub-set the data, add the argument \emph{subset.data=TRUE} and you 
#' can specify the number of cells to be used using \emph{n.samples} argument. For example \emph{n.samples=300}.
#' 
#' Different experimental designs can be simulated using the group and batch configuration arguments to
#' simulate biologica/experimental conditions and batchs (instrument or subjects), respectively. However,
#' the design of the new simulated data strongly depends on the design in the source data. In addition,
#' batch simulation may lead to unexpected error if the number of batches in the source data is small (e.g.
#' less than 5) and/or there is no enough data in each batch. 
#' 
#' 
#' @references
#' \itemize{
#' \item Efron, B., & Tibshirani, R. (1996). Using specially designed exponential families for density estimation. \emph{The Annals of Statistics}, 24(6), 2431-2461.
#' }
#' 
#' 
#' @examples 
#' ##
#' @export 
#' @importFrom plyr rbind.fill
simulateData <- function(n.sim=1, s.data, batch=NULL, group=NULL, n.genes=1000,  
                       batch.config= 1,  group.config=c(0.5, 0.5), 
                       tot.samples=150, pDE=0.2, model.zero.prob=FALSE, 
                       result.format="SCE", verbose=TRUE, ...)
{
  
  # experiment configurartion
  if(verbose) {message("Configuring design ...")}
  null.group = ifelse(!is.null(group), which.max(table(group))[[1]], 1)
  exprmt.design <- expriment.config(batch.config = batch.config, group.config = group.config,
                                    tot.samples = tot.samples)
  n.batch <- exprmt.design$n.batch
  n.group <- exprmt.design$n.group
  config.mat <- exprmt.design$exprmt.config 
  
  # prepare source data
  if(verbose) {message("Preparing source data ...")}
  prepare.S.Data <- prepareSourceData(s.data, batch = batch, group = group, 
                                      exprmt.design=exprmt.design,  ...)
  LL <- prepare.S.Data$LL
  cpm.data <- prepare.S.Data$cpm.data
  sub.batchs <- prepare.S.Data$sub.batchs
  
  if(is.null(group)) group <- rep(1, ncol(s.data))
  if(is.null(batch)) batch <- rep(1, ncol(s.data))
  
  # fit logistic regression for the probability of zeros
  if(model.zero.prob){
    if(verbose) {message("Fitting zero probability model ...")}
    if(class(s.data)=="SingleCellExperiment"){
      fracZero.logit.list <- lapply(unique(sub.batchs), function(b){
        zeroProbModel(counts(s.data)[, batch==b], ...)
      }) 
    }
    else if(class(s.data) %in% c("data.frame", "matrix")){
      fracZero.logit.list <- lapply(unique(sub.batchs), function(b){
        zeroProbModel(s.data[, batch==b], ...)
      })
    }
  }
  
  
  # candidate genes
  if(verbose) {message("Selecting genes ...")}
  null.genes0     <- prepare.S.Data$cand.genes$null.genes
  nonnull.genes0  <- prepare.S.Data$cand.genes$nonnull.genes
  
  
  
  # simulation step
  if(verbose) {message("Simulating data ...")}
  sim.data.list <- lapply(1:n.sim, function(h){
    if(verbose) {message(" ...", h, " of ", n.sim)}
    
    # sample DE and null genes
    if(pDE>0 & !is.null(group) & length(group.config)>1){
      if((1-pDE)*n.genes <= length(null.genes0)){
        null.genes     <- sample(null.genes0, (1-pDE)*(n.genes), replace = FALSE)
      }
      else{
        message("Note: The number of null genes (not DE) in the source data is ", length(null.genes0), 
                " and the number of null genes required to be included in the simulated data is ",
                round((1-pDE)*(n.genes)), ". Therefore, candidiate null genes are sampled with replacement.")
        null.genes     <- sample(null.genes0, (1-pDE)*(n.genes), replace = TRUE)
      }
      
      if(pDE*n.genes <= length(nonnull.genes0)){
        nonnull.genes  <- sample(nonnull.genes0, pDE*n.genes, replace = FALSE)
      }
      else{
        message("Note: The number of DE genes detected in the source data is ", length(nonnull.genes0), 
                " and the number of DE genes required to be included in the simulated data is ",
                round(pDE*n.genes), ". Therefore, candidiate DE genes are sampled with replacement.")
        nonnull.genes  <- sample(nonnull.genes0, pDE*n.genes, replace = TRUE)
      }
      
      sel.genes <- c(nonnull.genes, null.genes)
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
    
    # if(class(s.data)=="SingleCellExperiment"){
    #   sub.dat <- s.data[sel.genes, ]
    #   L <- colSums(counts(sub.dat))
    #   #cpm.data <- log2(calCPM(sub.dat)+1)
    # }
    # else if(class(s.data) %in% c("data.frame", "matrix")){
    #   sub.dat <- s.data[sel.genes, ]
    #   L <- colSums(sub.dat)
    #   #cpm.data <- log2(calCPM(sub.dat)+1)
    # }
    # LL <- lapply(1:length(sub.batchs), function(b){
    #   L.b <- L[batch==sub.batchs[b]]
    #   fit.ln <- fitdist(L.b, distr = "lnorm")$estimate
    #   L.b.pred <- rlnorm(n.batch[b], fit.ln[["meanlog"]], fit.ln[["sdlog"]])
    #   
    #   ## randomly split into the groups
    #   gr <- rep(1:length(n.group), config.mat[b, ])
    #   split(L.b.pred, gr)
    #   #LL.splited <- lapply(sort(unique(group)), function(g) L.b.pred[group==g & batch==b])
    #   #names(LL.splited) <- paste0("grp_", sort(unique(group)))
    #   #LL.splited
    # }) 
    
    
    # estimate batch specific parameters
    est.list <- lapply(sel.genes, function(i){
      #print(i)
      batch.est <- lapply(sort(unique(batch)), function(b){ 
        #print(b)
        if(DE.ind[i]==0){
          Y0 <- cpm.data[i, (batch==b & group==null.group)] 
          if(model.zero.prob & mean(Y0==0)>0.25) Y <- Y0[Y0>0]
          else Y <- Y0
          
          if(sum(Y>0)>3){
            countY   <- obtCount(Y)
            # plot(countY$S, countY$Ny, type="b") 
            parm.est <- fitLLmodel(countY)
            if(!is.null(parm.est$parm.list$betas) & length(countY$S)>=3){
              parm.est
            } 
            else{NULL}  
          } 
          else{NULL}
        }
        else{
          Y0 <- split(cpm.data[i, batch==b], group[batch==b])
          #Y0 <- lapply(sort(unique(group)), function(g) cpm.data[i, batch==b & group==g])
          #names(Y0) <-  paste0("grp_", sort(unique(group)))
          
          Y <- lapply(Y0, function(y0){
            if(model.zero.prob & mean(y0==0)>0.25) {y0[y0>0]}
            else {y0}
          }) 
          if(all(sapply(Y, function(y) sum(y>0)>3))){
            countY   <- lapply(Y,  obtCount)
            parm.est <- lapply(countY, fitLLmodel)
            
            cond1 = all(sapply(parm.est, function(x) !is.null(x$parm.list$betas)))
            cond2 = all(sapply(countY, function(x) length(x$S)>=3))
            if(cond1 & cond2){  
              parm.est
            } 
            else{NULL}
          } 
          else{NULL}
          # plot(countY$`1`$S, countY$`1`$Ny, type="l", xlim=range(cpm.data[i, batch==b])) 
          # lines(countY$`2`$S, countY$`2`$Ny, col=2)    
        }
      }) 
      
      if(DE.ind[i]==0){
        batch.parms.lst <- lapply(batch.est, function(b){ 
          if(!is.null(b)){
            c(b$parm.list$betas,  mu.hat=b$parm.list$mu.hat,  sig.hat=b$parm.list$sig.hat)
          }
          else{as.matrix(NA)}
        })
        batch.parms <- rbind.fill(lapply(batch.parms.lst, 
                                         function(x) as.data.frame(as.list(x)))) 
        if(any(dim(batch.parms)>1)){
          Mu.batch <- colMeans(batch.parms, na.rm = TRUE)
          V.batch  <- var(batch.parms, na.rm = TRUE) 
        }
        else{
          Mu.batch <- NULL
          V.batch  <- NULL
        }
      }
      else{
        batch.parms.lst <- lapply(sort(unique(group)), function(g){
          b <- lapply(batch.est, function(x) x[[g]])
          b <- lapply(b, function(bb){
            if(!is.null(bb)){
              c(bb$parm.list$betas, mu.hat=bb$parm.list$mu.hat, 
                sig.hat=bb$parm.list$sig.hat)
            }
            else{NULL}
          }) })
        
        batch.parms <- lapply(batch.parms.lst, function(g){
          rbind.fill(lapply(g, function(x) as.data.frame(as.list(x))))
        })
        
        if(!any(sapply(batch.parms, function(x) is.null(unlist(x))))){
          Mu.batch <- lapply(batch.parms, colMeans, na.rm = TRUE)
          V.batch  <- lapply(batch.parms, var, na.rm = TRUE)
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
      if(!is.null(est.list[[i]]) & !any(sapply(est.list[[i]]$batch.est, is.null))){
        # print(i)
        # sample params from MVN 
        if(DE.ind[i]==0){
          if(length(n.batch)>1){
            scl <- ifelse(n.batch<20, log10(log2(n.batch))+0.1, 1)
            par.sample <- mvrnorm(n = length(sub.batchs), mu= est.list[[i]]$Mu.batch,
                                  Sigma = scl*est.list[[i]]$V.batch) 
          }
          else if(length(n.batch)==1){
            par.sample <- t(as.matrix(est.list[[i]]$Mu.batch))
          } 
        }
        else{ 
          if(length(n.batch)>1){
            scl <- ifelse(n.batch<20, log10(log2(n.batch))+0.1, 1)
            par.sample <- lapply(sort(unique(group)), function(g){
              mvrnorm(n = length(sub.batchs), mu= est.list[[i]]$Mu.batch[[g]],
                      Sigma = scl*est.list[[i]]$V.batch[[g]]) 
            })
            #names(par.sample) <- paste0("grp_", sort(unique(group)))
          }
          else if(length(n.batch)==1){
            par.sample <- lapply(sort(unique(group)), function(g){
              t(as.matrix(est.list[[i]]$Mu.batch[[g]]))
            })
            #names(par.sample) <- paste0("grp_", sort(unique(group)))
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
          g0 <- lapply(sort(unique(group)), function(g){
            par.sample.g <- par.sample[[g]] #par.sample[[paste0("grp_", g)]]
            lapply(1:nrow(par.sample.g), function(bb){  
              b.data <- est.list[[i]]$batch.est[[bb]][[g]]
              gg0 <- (pnorm(b.data$yy$uls, par.sample.g[bb, "mu.hat"][[1]], par.sample.g[bb, "sig.hat"][[1]]) -
                        pnorm(b.data$yy$lls, par.sample.g[bb, "mu.hat"][[1]], par.sample.g[bb, "sig.hat"][[1]]))
              gg0[is.nan(gg0)] <- 0
              gg0 
            })
          })
          #names(g0) <- paste0("grp_", sort(unique(group)))
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
          g1 <- lapply(sort(unique(group)), function(g){
            par.sample.g <- par.sample[[g]] #par.sample[[paste0("grp_", g)]]
            g0.g <- g0[[g]] #g0[[paste0("grp_", g)]]
            
            lapply(1:nrow(par.sample.g), function(bb){   
              b.data <- est.list[[i]]$batch.est[[bb]][[g]]
              gg0 <- g0.g[[bb]]*sum(b.data$yy$Ny) + 1
              
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
          #names(g1) <- paste0("grp_", sort(unique(group)))
        }
        
        # simulate data 
        if(DE.ind[i]==0){ 
          Y.star <- lapply(1:nrow(par.sample), function(bb){
            Y0 <- cpm.data[sel.genes[i], (batch==bb & group==null.group)] 
            u <- runif(n.batch[[bb]])
            gg1 <- g1[[bb]]
            y.star.b <- sapply(u, function(uu){
              yy  <- gg1$s[which.min(abs(gg1$Gy-uu))]
              eps <- gg1$s[2]-gg1$s[1] #=mean(diff(gg1$s))
              yy  <- runif(1, yy-eps/2, yy+eps/2)
              yy
            })
            LL.b <- as.numeric(do.call("c", LL[[bb]]))
            y.star.b <- round(((exp(y.star.b)-1)*LL.b)/1e6)
            y.star.b[y.star.b<0] <- 0
            
            if(model.zero.prob & mean(Y0==0)>0.25){
              lLL_b <- log(LL.b)
              pred.pz  <- predict(fracZero.logit.list[[bb]], type="response",
                                  newdata=data.frame(x1=mean(Y0), x2=lLL_b))
              drop.mlt <- sapply(pred.pz, function(p) rbinom(1, 1, p))
              y.star.b <- y.star.b*(1-drop.mlt)
            }
            as.numeric(y.star.b)
          })
          Y.star <- do.call("c", Y.star)
          #if(any(is.na(Y.star))){print(i)}
        }
        else{ 
          Y.star <- lapply(sort(unique(group)), function(g){
            par.sample.g <- par.sample[[g]] #par.sample[[paste0("grp_", g)]]
            g1.g <- g1[[g]] #g1[[paste0("grp_", g)]]
            y.star.g <- lapply(1:nrow(par.sample.g), function(bb){
              Y0.g <- cpm.data[sel.genes[i], batch==bb & group==g]
              u <- runif(config.mat[bb, g])#runif(n.batch[[bb]]/length(n.group))
              gg1 <- g1.g[[bb]]
              y.star.b <- sapply(u, function(uu){
                yy  <- gg1$s[which.min(abs(gg1$Gy-uu))]
                eps <- gg1$s[2]-gg1$s[1] #=mean(diff(gg1$s))
                yy  <- runif(1, yy-eps/2, yy+eps/2)
                yy
              })
              LL.b.g <- LL[[bb]][[g]]
              y.star.b <- round(((exp(y.star.b)-1)*LL.b.g)/1e6)
              y.star.b[y.star.b<0] <- 0
              
              if(model.zero.prob & mean(Y0.g==0)>0.25){
                lLL.b.g<- log(LL.b.g)
                pred.pz  <- predict(fracZero.logit.list[[bb]], type="response",
                                    newdata=data.frame(x1=mean(Y0.g), x2=lLL.b.g))
                drop.mlt <- sapply(pred.pz, function(p) rbinom(1, 1, p))
                y.star.b <- y.star.b*(1-drop.mlt)
              }
              as.numeric(y.star.b)
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
                           sim.Lib.Size=do.call("c", do.call("c", LL)),
                           LS = colSums(sim.count, na.rm=TRUE), 
                           row.names = colnames(sim.count))
    row.data <- data.frame(DE.ind=DE.ind, 
                           source.ID=sel.genes,
                           row.names = rownames(sim.count))
    
    if(result.format == "SCE"){
      sim.data <- SingleCellExperiment(assays=list(counts=sim.count),
                                       colData=col.data, rowData=row.data)
      sim.data 
    }
    else if(result.format == "list"){
      sim.data <- list(counts=sim.count, colData=col.data, rowData=row.data)
      sim.data
    }
  })
  sim.data.list
}

