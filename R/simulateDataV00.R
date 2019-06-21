# A function to simulate bulk or single cell RNA sequencing data
# 
# @description This function simulates RNA sequencing data given a real RNA-seq data using
# semi-parametric density estimation.
# 
# @param n.sim a numerical value for the number of simulated data to be  generated
# @param s.data a source data (a SingleCellExperiment class object or a matrix/data.frame of counts with genes in 
# rows and samples in columns)
# @param batch a vector containing btach indicator for each sample/cell
# @param group a vector containg group indicator for each sample/cell 
# @param n.genes a numeric value for the total number of genes to be simulated
# @param pDE a numeric value between 0 and 1 indicating the fraction of DE genes 
# in a single simulated data
# @param batch.config a numerical vector for the marginal fraction of samples in each batch. 
# The number of batches to be simulated is equal to the size of the vector.
# All values must sum to 1.
# @param group.config a numerical vector for the marginal fraction of samples in each group.
# The number of groups to be simulated is equal to the size of the vector. All values must sum to 1.
# @param model.zero.prob a logical value whether to model the zero probablity separately 
# (suitable for single cell data)
# @param tot.samples a numerical value for total number of samples to be simulated. 
# @param result.format a character value for the type of format for the output. Choice can  be
# 'SCE' for SingleCellExperiment class or "list" for a list object that contains the simulated count,
# column information abd row information.
# @param const a small constant (>0) to be added to the CPM before log transformation, to avoid  log(0).
# default is 1e-5
# @param verbose a logical value, if TRUE it displays a message about the satatus of the simulation
# @param  seed an integer  between 1 and 1e10. It will be used for #set.seed() function
# @param  ... further arguments passed to or from other methods.
# 
# @return a list of SingleCellExperiment object each contatining simulated counts (not normalized), 
# cell level information in colData, and gene level information in rowData.
# 
# 
# @details This function estimates the density of a given bulk or single cell RNA-seq data 
# (passed using \emph{s.data} argument) using a specially designed exponetial family for density
# estimation. Afterwards, it simulates a new data set from the estimated density. In a first step, 
# the log-CPM outcomes from a given real data  are used for semi-parametrically 
# estimating gene-wise distributions. This method is based on a fast log-linear model estimation 
# approach developed by Efron et al (1996). Arbitrarily large datasets, with realistically varying 
# library sizes, can be sampled from these estimated distributions.
# 
# For simulation of single cell RNA-seq data, it involves an additional step to explicitly account for the high 
# abundance of zero counts. This step models the probability of zero counts as a function the mean 
# expression of the gene and the library size of the cell (both in log scale) to add excess zeros. 
# This can be done by using \emph{model.zero.prob=TRUE}. Note that,
# for extremly large size data, it is recomended to use a random sample of cells to
# reduce computation time. To enable this, add the argument \emph{subset.data=TRUE} and you 
# can specify the number of cells to be used using \emph{n.samples} argument. 
# For example \emph{n.samples=400}.
# 
# Given known groups of samples/cells in the source data, DGE is simulated by independently 
# sampling data from distributions constructed in each group. In particular, this procedure is 
# applied on a set of genes with fold-change in the source data more than a given threshold (\emph{lfc.thrld}). 
# Moreover, when  the source dataset involves samples/cells processed in different batches, our 
# simulation procedure incorporates this batch effect in the simulated data, if required. 
# Different experimental designs can be simulated using the group and batch configuration arguments to
# simulate biologica/experimental conditions and batchs (instrument or subjects), respectively. 
# Also, it is important to filter the source data so that gene with suffient expression will be used to 
# estimate the density.
# 
# 
# @references
# \itemize{
# \item Efron, B., & Tibshirani, R. (1996). Using specially designed exponential families for density estimation. \emph{The Annals of Statistics}, 24(6), 2431-2461.
# }
# 
# 
# @examples 
# #----------------------------------------------------------------
# # Example 1: simulating bulk RNA-seq
# 
# \donttest{ # TAKES LONG
# # load the Zhang data (availabl with the package)
# data("zhang.data") 
# 
# # filter genes with sufficient expression (important step to avoid bugs) 
# zhang.counts <- zhang.data$counts[rowSums(zhang.data$counts > 0)>=10, ]  
# MYCN.status <- zhang.data$MYCN.status+1  
# 
# # We simulate only a single data (n.sim = 1) with the following property
# # - 2000 genes ( n.genes = 2000) 
# # - 180 samples (tot.samples = 180) 
# # - the samples are equally divided into 2 groups each with 90 samples 
# #   (group.config = c(0.5, 0.5))
# # - all samples are from a single batch (batch.config = 1)
# # - we add 10% DE genes (pDE = 0.1) 
# # - we do not model the zeroes separately, they are the part of density 
# #    estimation (model.zero.prob = FALSE)
# 
# sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = zhang.counts, batch = NULL,
#                               group = MYCN.status, n.genes = 2000, batch.config = 1,
#                               group.config = c(0.5, 0.5), tot.samples = 180, pDE = 0.1,
#                               model.zero.prob = FALSE, result.format = "list")
#                               
# sim.data.bulk1 <- sim.data.bulk[[1]]                              
# head(sim.data.bulk1$counts[, 1:5])  # count data
# head(sim.data.bulk1$colData)        # sample info
# head(sim.data.bulk1$rowData)        # gene info
# 
# 
# #----------------------------------------------------------------
# # Example 2: simulating single cell RNA-seq from a single batch (read-counts)
# # we simulate only a single scRNA-seq data (n.sim = 1) with the following property
# # - 2000 genes (n.genes = 2000) 
# # - 100 cells (tot.samples = 100) 
# # - the cells are equally divided into 2 groups each with 50 cells (group.config = c(0.5, 0.5))
# # - all cells are from a single batch (batch.config = 1)
# # - we add 10% DE genes (pDE = 0.1) 
# # - we model the zeroes separately (model.zero.prob = TRUE)
# # - the ouput will be in SingleCellExperiment class object (result.format = "SCE")
# 
# 
# library(SingleCellExperiment)
# 
# # load the NGP nutlin data (availabl with the package)
# data("scNGP.data")
# 
# # filter genes with sufficient expression (important step to avoid bugs) 
# scNGP.data2 <- scNGP.data[rowSums(counts(scNGP.data) > 0)>=10, ]  
# treatment <- ifelse(scNGP.data2$characteristics..treatment=="nutlin",2,1) 
# 
# # simulate data (we simulate here only a single data, n.sim = 1)
# sim.data.sc <- SPsimSeq(n.sim = 1, s.data = scNGP.data2, batch = NULL,
#                             group = treatment, n.genes = 2000, batch.config = 1,
#                             group.config = c(0.5, 0.5), tot.samples = 100, pDE = 0.1,
#                             model.zero.prob = TRUE, result.format = "SCE")
#                             
# sim.data.sc1 <- sim.data.sc[[1]]
# class(sim.data.sc1)
# head(counts(sim.data.sc1)[, 1:5])
# colData(sim.data.sc1)
# rowData(sim.data.sc1)
# 
# 
# 
# #----------------------------------------------------------------
# # Example 3: simulating single cell RNA-seq from a single batch (UMI counts)
# # we simulate only a single scRNA-seq data (n.sim = 1) with the following property
# # - 2000 genes (n.genes = 2000) 
# # - 200 cells (tot.samples = 200) 
# # - the cells are from a single experimental group (group.config = 1)
# # - all cells are from a single batch (batch.config = 1)
# # - we add 0% DE genes (pDE = 0) 
# # - we model the zeroes separately (model.zero.prob = TRUE)
# # - since the size of the PBMC data is large, we use the subset of the cells to 
# #   fit the zero prob. model (subset.data=TRUE, n.samples=400)
# # - the ouput will be in SingleCellExperiment class object (result.format = "SCE") 
# 
# library(SingleCellExperiment)
# 
# # load the PBMC data (availabl with the package)
# data("PBMC.data") 
# 
# # filter genes with sufficient expression (important step to avoid bugs) 
# PBMCdat2 <- PBMC.10x.data[rowSums(counts(PBMC.10x.data) > 0)>=20, ] 
# 
# # simulate data (we simulate here only a single data, n.sim = 1)
# sim.data.scUMI <- SPsimSeq(n.sim = 1, s.data = PBMCdat2, batch = NULL,
#                                group = NULL, n.genes = 2000, batch.config = 1,
#                                group.config = 1, tot.samples = 200, pDE = 0,
#                                model.zero.prob = TRUE, result.format = "SCE", 
#                                subset.data=TRUE, n.samples=400)
#                             
# sim.data.scUMI1 <- sim.data.scUMI[[1]]
# class(sim.data.scUMI1)
# head(counts(sim.data.scUMI1)[, 1:5])
# colData(sim.data.scUMI1)
# rowData(sim.data.scUMI1)
# 
# }
# 
# 
# @export    
# @importFrom MASS mvrnorm
# @importFrom stats pnorm runif rbinom predict approx quantile glm rlnorm rnbinom 
# @importFrom SingleCellExperiment counts colData rowData SingleCellExperiment
# @importFrom Hmisc cut2
# @importFrom fitdistrplus fitdist
# @importFrom stats glm pnorm coef vcov
# @importFrom edgeR calcNormFactors
# @importFrom SingleCellExperiment counts SingleCellExperiment colData rowData
# @importFrom graphics hist
.SPsimSeqV00 <- function(n.sim=1, s.data, batch=NULL, group=NULL, n.genes=1000,  
                       batch.config= 1,  group.config=c(0.5, 0.5), 
                       tot.samples=150, pDE=0.2, model.zero.prob=FALSE, const=1e-5,
                       result.format="SCE", verbose=TRUE,  seed=2581988, ...)
{
  # experiment configurartion
  if(verbose) {message("Configuring design ...")}
  null.group = ifelse(!is.null(group), which.max(table(group))[[1]], 1)
  exprmt.design <- expriment.config(batch.config = batch.config, group.config = group.config,
                                    tot.samples = tot.samples)
  n.batch <- exprmt.design$n.batch
  n.group <- exprmt.design$n.group
  config.mat <- exprmt.design$exprmt.config 
  
  # sim control
  #setSimContol <- setSimContol(seed, n.sim = n.sim, n.genes = n.genes)
  if(!is.null(seed)){set.seed(seed)} 
  
  # prepare source data
  if(verbose) {message("Preparing source data ...")}
  prepare.S.Data <- prepareSourceData(s.data, batch = batch, group = group, 
                                      exprmt.design=exprmt.design, const=const,
                                      simCtr=NULL, ...)
  LL         <- prepare.S.Data$LL
  cpm.data   <- prepare.S.Data$cpm.data
  sub.batchs <- prepare.S.Data$sub.batchs
  
  if(is.null(group)) group <- rep(1, ncol(s.data))
  if(is.null(batch)) batch <- rep(1, ncol(s.data))
  
  # fit logistic regression for the probability of zeros
  if(model.zero.prob){
    if(verbose) {message("Fitting zero probability model ...")}
    if(class(s.data)=="SingleCellExperiment"){
      fracZero.logit.list <- lapply(unique(sub.batchs), function(b){
        zeroProbModel(cpm.data = cpm.data[, batch==b], L=colSums(counts(s.data)[, batch==b]),
                      simCtr=NULL, const=const, ...)
      }) 
    }
    else if(class(s.data) %in% c("data.frame", "matrix")){
      fracZero.logit.list <- lapply(unique(sub.batchs), function(b){
        zeroProbModel(cpm.data = cpm.data[, batch==b], L=colSums(s.data[, batch==b]), 
                      simCtr=NULL, const,...)
      })
    }
  }
  else{
    fracZero.logit.list=NULL
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
    selctGenes <- selectGenesSim(pDE = pDE, group = group, n.genes = n.genes, null.genes0 = null.genes0,
                  nonnull.genes0 = nonnull.genes0, group.config = group.config)
    DE.ind <- selctGenes$DE.ind
    sel.genes <- selctGenes$sel.genes
    
    
    
    # estimate batch specific parameters
    min.val <- log(const)
    est.list <- lapply(sel.genes, function(i){
      #print(i)
      gene.parm.est(cpm.data.i = cpm.data[i, ], batch = batch, group = group, 
                    null.group = null.group, sub.batchs = sub.batchs, de.ind = DE.ind[i], 
                    model.zero.prob = model.zero.prob, min.val = min.val, ...)
    })
    
    sim.list <- lapply(1:length(sel.genes), function(i){
      #print(i)
      if(!is.null(est.list[[i]]) & 
         !any(sapply(est.list[[i]]$batch.est, is.null))) # & 
         #!any(is.na(est.list[[i]]$V.batch)))
        {
        # print(i)
        # sample params from MVN 
        if(DE.ind[i]==0){
          if(length(n.batch)>1){  
             
            #parm.seed <- setSimContol$seed.batch.pars[i]
            par.sample <- as.matrix(do.call("rbind.fill2", lapply(est.list[[i]]$batch.est, function(bt){
              v.mat     <- bt$parm.list$v
              betas.vec <- bt$parm.list$betas 
              data.frame(t(as.matrix(c(mvrnorm(n = 1, mu= betas.vec, Sigma = v.mat), 
                mu.hat=bt$parm.list$mu.hat, sig.hat=bt$parm.list$sig.hat))))
            }))) 
          }
          else if(length(n.batch)==1){
            par.sample <- t(as.matrix(est.list[[i]]$Mu.batch)) 
          } 
        }
        else{ 
          if(length(n.batch)>1){
            #parm.seed <- setSimContol$seed.batch.pars[i]
            par.sample <- lapply(sort(unique(group)), function(g){ 
              par.sample.g <- as.matrix(do.call("rbind.fill2", 
                            lapply(est.list[[i]]$batch.est[[g]], function(bt){
                              #set.seed(parm.seed)
                              data.frame(t(as.matrix(c(mvrnorm(n = 1, mu= bt$parm.list$betas, 
                                                               Sigma = bt$parm.list$v), 
                                                       mu.hat=bt$parm.list$mu.hat,
                                                       sig.hat=bt$parm.list$sig.hat)))) 
              })))
              
              par.sample.g 
            })
            #set.seed(NULL)
            #names(par.sample) <- paste0("grp_", sort(unique(group)))
          }
          else if(length(n.batch)==1){
            par.sample <- lapply(sort(unique(group)), function(g){
              par.sample.g <-t(as.matrix(est.list[[i]]$Mu.batch[[g]])) 
              par.sample.g 
            }) 
          } 
        }
        
        
        # estimate carrier density (g0)
        g0 <- SPestCarrier(de.ind = DE.ind[i], par.sample = par.sample, group=group, 
                           est.list.i = est.list[[i]], ...)
        # if(DE.ind[i]==0){
        #   g0 <- lapply(1:nrow(par.sample), function(bb){
        #     b.data <- est.list[[i]]$batch.est[[bb]] 
        #     gg0 <- (pnorm(b.data$yy$uls, as.matrix(par.sample)[bb, "mu.hat"],
        #                   as.matrix(par.sample)[bb, "sig.hat"]) -
        #               pnorm(b.data$yy$lls,as.matrix(par.sample)[bb, "mu.hat"], 
        #                     as.matrix(par.sample)[bb, "sig.hat"]))
        #     
        #     gg0[is.nan(gg0)] <- 0
        #     gg0.expd <- sapply(b.data$yy$Y, function(y) gg0[b.data$yy$lls<y & b.data$yy$uls>=y]) 
        #     gg0.expd
        #   })
        # }
        # else{
        #   g0 <- lapply(sort(unique(group)), function(g){
        #     par.sample.g <- par.sample[[g]] #par.sample[[paste0("grp_", g)]]
        #     lapply(1:nrow(par.sample.g), function(bb){  
        #       b.data <- est.list[[i]]$batch.est[[bb]][[g]] 
        #       gg0 <- (pnorm(b.data$yy$uls, par.sample.g[bb, "mu.hat"][[1]], 
        #                     par.sample.g[bb, "sig.hat"][[1]]) -
        #                 pnorm(b.data$yy$lls, par.sample.g[bb, "mu.hat"][[1]], 
        #                       par.sample.g[bb, "sig.hat"][[1]]))
        #       gg0[is.nan(gg0)] <- 0 
        #       gg0.expd <- sapply(b.data$yy$Y, function(y) gg0[b.data$yy$lls<y & b.data$yy$uls>=y]) 
        #       gg0.expd
        #     })
        #   })
        #   #names(g0) <- paste0("grp_", sort(unique(group)))
        # }
        
        
        
        # estimate the density g1(y)
        g1 <- SPestDens(de.ind = DE.ind[i], group = group, batch = batch, 
                        par.sample = par.sample, est.list.i = est.list[[i]], g0 = g0, ...)
        # if(DE.ind[i]==0){
        #   g1 <- lapply(1:nrow(par.sample), function(bb){ 
        #     b.data <- est.list[[i]]$batch.est[[bb]]
        #     gg0 <- g0[[bb]]*sum(b.data$yy$Ny) +1
        #     SPestDens(b.data = b.data, gg0 = gg0, par.sample = par.sample, bb = bb)
        #     
        #     # s <- b.data$yy$Y
        #     # s.mat <- matrix(NA, ncol = length(coef(b.data$llm)), nrow=length(s))
        #     # for(ii in 1:ncol(s.mat)){
        #     #   s.mat[, ii] <- s^(ii-1)
        #     # }
        #     # gg1 <- exp(s.mat %*% par.sample[bb, 1:ncol(s.mat)])*gg0
        #     # #gg1 <- data.frame(gy=gg1, s=s, lls=b.data$yy$lls, uls=b.data$yy$uls)
        #     # lls.expd <- sapply(s, function(ss) b.data$yy$lls[b.data$yy$lls<ss & b.data$yy$uls>=ss])
        #     # uls.expd <- sapply(s, function(ss) b.data$yy$uls[b.data$yy$lls<ss & b.data$yy$uls>=ss])
        #     # S.expd   <- sapply(s, function(ss) b.data$yy$S[b.data$yy$lls<ss & b.data$yy$uls>=ss])
        #     # Ny.expd  <- sapply(s, function(ss) b.data$yy$Ny[b.data$yy$lls<ss & b.data$yy$uls>=ss])
        #     # gg1 <- data.frame(s=s, S.expd=S.expd, lls.expd=lls.expd, 
        #     #                   uls.expd=uls.expd, Ny.expd=Ny.expd, gg0=gg0, gy=gg1)
        #     # gg1 <- gg1[!duplicated(gg1),]
        #     # gg1$gy[is.infinite(gg1$gy)] <- max(gg1$gy[!is.infinite(gg1$gy)]) 
        #     # gg1 <- gg1[order(gg1$s),]
        #     # gg1$Gy <- cumsum(gg1$gy/gg1$Ny.expd)/sum(gg1$gy/gg1$Ny.expd) 
        #     # gg1
        #     
        #     
        #   })  
        #   # plot(g1[[1]]$s, g1[[1]]$gy, type="b", xlim=range(do.call('c', lapply(g1, function(x) x[, "s"]))))
        #   # for(k in 2:length(g1)){
        #   #   lines(g1[[k]]$s, g1[[k]]$gy, type="b", col=k)
        #   # }
        # }
        # else{
        #   g1 <- lapply(sort(unique(group)), function(g){
        #     par.sample.g <- par.sample[[g]] #par.sample[[paste0("grp_", g)]]
        #     g0.g <- g0[[g]] #g0[[paste0("grp_", g)]]
        #     
        #     lapply(1:nrow(par.sample.g), function(bb){   
        #       b.data <- est.list[[i]]$batch.est[[bb]][[g]]
        #       gg0 <- g0.g[[bb]]*sum(b.data$yy$Ny) + 1
        #       
        #       SPestDens(b.data = b.data, gg0 = gg0, par.sample = par.sample, bb = bb)
        #       # s <- b.data$yy$Y 
        #       # s.mat <- matrix(NA, ncol = length(coef(b.data$llm)), nrow=length(s))
        #       # for(k in 1:ncol(s.mat)){
        #       #   s.mat[, k] <- s^(k-1)
        #       # }
        #       # gg1 <- exp(s.mat %*% as.matrix(par.sample.g[bb, 1:ncol(s.mat)]))*gg0 
        #       # #gg1 <- data.frame(gy=gg1, s=s, lls=b.data$yy$lls, uls=b.data$yy$uls)
        #       # lls.expd <- sapply(s, function(ss) b.data$yy$lls[b.data$yy$lls<ss & b.data$yy$uls>=ss])
        #       # uls.expd <- sapply(s, function(ss) b.data$yy$uls[b.data$yy$lls<ss & b.data$yy$uls>=ss])
        #       # S.expd   <- sapply(s, function(ss) b.data$yy$S[b.data$yy$lls<ss & b.data$yy$uls>=ss])
        #       # Ny.expd  <- sapply(s, function(ss) b.data$yy$Ny[b.data$yy$lls<ss & b.data$yy$uls>=ss])
        #       # 
        #       # gg1 <- data.frame(s=s, S.expd=S.expd, lls.expd=lls.expd, 
        #       #                   uls.expd=uls.expd, Ny.expd=Ny.expd, gg0=gg0, gy=gg1)
        #       # gg1 <- gg1[!duplicated(gg1),]
        #       # gg1$gy[is.infinite(gg1$gy)] <- max(gg1$gy[!is.infinite(gg1$gy)]) 
        #       # gg1 <- gg1[order(gg1$s),]
        #       # gg1$Gy <- cumsum(gg1$gy/gg1$Ny.expd)/sum(gg1$gy/gg1$Ny.expd) 
        #       # gg1
        #     }) 
        #   })
        #   #names(g1) <- paste0("grp_", sort(unique(group)))
        # }
        
        # simulate data 
        #sim.seed <- setSimContol$seed.sample.Y[i]
        Y.star <- SPsampleData(de.ind = DE.ind[i], par.sample = par.sample, cpm.data.i = cpm.data[i,], 
                               batch =batch, group = group, null.group = null.group, g1 = g1, LL = LL,
                               model.zero.prob = model.zero.prob, min.val = min.val, const = const,
                               config.mat = config.mat, n.batch=n.batch, 
                               fracZero.logit.list=fracZero.logit.list, ...)
        # if(DE.ind[i]==0){ 
        #   Y.star <- lapply(1:nrow(par.sample), function(bb){
        #     Y0 <- cpm.data[sel.genes[i], (batch==bb & group==null.group)] 
        #     #set.seed(sim.seed)
        #     u <- runif(n.batch[[bb]]) 
        #     gg1 <- g1[[bb]]
        #     rownames(gg1) <- 1:nrow(gg1)
        #     #set.seed(sim.seed+1)
        #     y.star.b <- sapply(u, function(uu){
        #       whc.r  <- sample(which.min(abs(gg1$Gy-uu)), 1) # incase there are multiples
        #       yy     <- gg1$s[whc.r] 
        #       # difs<- diff(gg1$s)
        #       # difs<- difs[!(is.na(difs) | is.infinite(difs) | is.nan(difs))]
        #       # eps <- abs(mean(difs, na.rm=TRUE)) #gg1$s[2]-gg1$s[1] 
        #       if(whc.r==1){
        #         yll <- max(yy - 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #                    gg1$lls.expd[whc.r])
        #         yul <- min(min(yy + 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #                        gg1$uls.expd[whc.r]), gg1$s[whc.r+1])
        #       }else if(whc.r==nrow(gg1)){
        #         yll <- max(gg1$s[whc.r-1], max(yy - 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #                                        gg1$lls.expd[whc.r]))
        #         yul <- min(yy + 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #                    gg1$uls.expd[whc.r])
        #       }else{
        #         yll <- max(gg1$s[whc.r-1], max(yy - 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #                                        gg1$lls.expd[whc.r]))
        #         yul <- min(min(yy + 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #                        gg1$uls.expd[whc.r]), gg1$s[whc.r+1])
        #       }
        #       yy  <- suppressWarnings(runif(1, yll, yul))
        #       yy
        #     })
        #     LL.b <- as.numeric(do.call("c", LL[[bb]]))
        #     y.star.b <- round(((exp(y.star.b)-const)*LL.b)/1e6)
        #     y.star.b[y.star.b<0] <- 0
        #     
        #     if(model.zero.prob & mean(Y0==min.val)>0.25){
        #       lLL_b <- log(LL.b)
        #       pred.pz  <- try(predict(fracZero.logit.list[[bb]], type="response",
        #                           newdata=data.frame(x1=mean(Y0), x2=lLL_b)), 
        #                       silent = TRUE)
        #       if(class(pred.pz) != "try-error"){
        #         #set.seed(sim.seed+2)
        #         drop.mlt <- sapply(pred.pz, function(p){ 
        #           rbinom(1, 1, p) 
        #           })
        #       }
        #       else{
        #         drop.mlt <- 0
        #       } 
        #       y.star.b <- y.star.b*(1-drop.mlt)
        #     }
        #     as.numeric(y.star.b)
        #   })
        #   Y.star <- do.call("c", Y.star)
        #   #if(any(is.na(Y.star))){print(i)}
        # }
        # else{ 
        #   Y.star <- lapply(sort(unique(group)), function(g){
        #     par.sample.g <- par.sample[[g]] #par.sample[[paste0("grp_", g)]]
        #     g1.g <- g1[[g]] #g1[[paste0("grp_", g)]]
        #     y.star.g <- lapply(1:nrow(par.sample.g), function(bb){
        #       Y0.g <- cpm.data[sel.genes[i], batch==bb & group==g]
        #       #set.seed(sim.seed+bb)
        #       u <- runif(config.mat[bb, g])#runif(n.batch[[bb]]/length(n.group))
        #       gg1 <- g1.g[[bb]]
        #       rownames(gg1) <- 1:nrow(gg1)
        #       #set.seed(sim.seed+bb+1)
        #       y.star.b <- sapply(u, function(uu){
        #         whc.r  <- sample(which.min(abs(gg1$Gy-uu)), 1) # incase there are multiples
        #         yy     <- gg1$s[whc.r]
        #         #difs<- diff(gg1$s)
        #         #difs<- difs[!(is.na(difs) | is.infinite(difs) | is.nan(difs))]
        #         #eps <- abs(mean(difs, na.rm=TRUE)) #gg1$s[2]-gg1$s[1]   
        #         if(whc.r==1){
        #           yll <- max(yy - 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #                      gg1$lls.expd[whc.r])
        #           yul <- min(min(yy + 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #                      gg1$uls.expd[whc.r]), gg1$s[whc.r+1])
        #         }else if(whc.r==nrow(gg1)){
        #           yll <- max(gg1$s[whc.r-1], max(yy - 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #                      gg1$lls.expd[whc.r]))
        #           yul <- min(yy + 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #                          gg1$uls.expd[whc.r])
        #         }else{
        #           yll <- max(gg1$s[whc.r-1], max(yy - 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #                                          gg1$lls.expd[whc.r]))
        #           yul <- min(min(yy + 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #                          gg1$uls.expd[whc.r]), gg1$s[whc.r+1])
        #         }
        #         yy  <- suppressWarnings(runif(1, yll, yul))
        #         yy
        #       })
        #       LL.b.g <- LL[[bb]][[g]]
        #       y.star.b <- round(((exp(y.star.b)-const)*LL.b.g)/1e6)
        #       y.star.b[y.star.b<0] <- 0
        #       
        #       if(model.zero.prob & mean(Y0.g==min.val)>0.25){
        #         lLL.b.g<- log(LL.b.g)
        #         pred.pz  <- try(predict(fracZero.logit.list[[bb]], type="response",
        #                             newdata=data.frame(x1=mean(Y0.g), x2=lLL.b.g)), 
        #                         silent = TRUE)
        #         if(class(pred.pz) != "try-error"){
        #           #set.seed(sim.seed+bb+2)
        #           drop.mlt <- sapply(pred.pz, function(p){
        #             ##set.seed(sim.seed+bb+2)
        #             rbinom(1, 1, p)
        #             })
        #         }
        #         else{
        #           drop.mlt <- 0
        #         } 
        #         y.star.b <- y.star.b*(1-drop.mlt)
        #       }
        #       as.numeric(y.star.b)
        #     })
        #     y.star.g <- do.call("c", y.star.g)
        #   })
        #   Y.star <- do.call("c", Y.star)
        #   #if(any(is.na(Y.star))){print(i)}
        # }
        }
      else{
        #print(i)
        Y.star <- rep(0, tot.samples)
        Y.star
      }
    }) 
    
    sim.count <- do.call(rbind, sim.list)
    sim.count[is.na(sim.count)] <- 0
    rownames(sim.count) <- paste0("Gene_", 1:nrow(sim.count))
    colnames(sim.count) <- paste0("Sample_", 1:ncol(sim.count))
     
    col.data <- data.frame(Batch=rep(rep(1:length(n.batch), times=length(n.group)), config.mat),
                           Group=rep(rep(1:length(n.group), each=length(n.batch)), config.mat), 
                           sim.Lib.Size=do.call("c", do.call("c", LL)),  
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

