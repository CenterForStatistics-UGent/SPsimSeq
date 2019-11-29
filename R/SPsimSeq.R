#' A function to simulate bulk or single cell RNA sequencing data
#' 
#' @description This function simulates (bulk/single cell) RNA-seq dataset from semi-parametrically estimated distributions of gene expression levels in a given real source RNA-seq dataset
#'  
#' @param n.sim an integer for the number of simulations to be generated
#' @param s.data a real source dataset (a SingleCellExperiment class or a matrix/data.frame of counts with genes in rows and samples in columns)
#' @param batch NULL or a vector containing batch indicators for each sample/cell in the source data
#' @param group NULL or a vector containing group indicators for each sample/cell in the source data
#' @param n.genes a numeric value for the total number of genes to be simulated
#' @param pDE a numeric value between 0 and 1 indicating the desired fraction of DE genes in the simulated data
#' @param batch.config a numerical vector containing fractions  for the composition of samples/cells per batch. The fractions must sum to 1. The number of batches to be simulated is equal to the size of the vector. (Example, batch.config=c(0.6, 0.4)   means simulate 2 batches with 60\% of the simulated samples/cells in batch 1 and the rest 40\% in the second batch. Another example,  batch.config=c(0.3, 0.35, 0.25) means simulate 3 batches with the first, second and third batches contain 30\%, 35\% and 25\% of the samples/cells, respectively). 
#' @param group.config a numerical vector containing fractions  for the composition of samples/cells per group. Similar definition to 'batch.config'. The number of groups to be simulated is equal to the size of the vector. The fractions must sum to 1.
#' @param lfc.thrld a positive numeric value for the minimum absolute log-fold-change for selecting candidate DE genes in the source data (when group is not NULL and pDE>0)
#' @param t.thrld a positive numeric value for the minimum absolute t-test statistic for the log-fold-changes of genes for selecting candidate DE genes in the source data (when group is not NULL and  pDE>0)
#' @param llStat.thrld a positive numeric value for the minimum squared test statistics from the log-linear model to select candidate DE genes in the source data (when group is not NULL and  pDE>0)
#' @param model.zero.prob a logical value whether to model the zero expression probability separately (suitable for simulating of single-cell RNA-seq data or zero-inflated data)
#' @param tot.samples a numerical value for the total number of samples/cells to be simulated. 
#' @param result.format a character value for the type of format for the output. Choice can be 'SCE' for SingleCellExperiment class or "list" for a list object that contains the simulated count, column information and row information.
#' @param geneTogenesCor a logical value, if TRUE (default) the simulation will retain the gene-to-gene correlation structure of the source data using Gausian-copulas . Note that if it is TRUE, the program will be slow or it may fail for a limited memory size. 
#' @param w a numeric value between 0 and 1. The number of classes to construct the probability distribution will be round(w*n), where n is the total number of samples/cells in a particular batch of the source data
#' @param const a positive constant to be added to the CPM before log transformation, to avoid log(0). The default is 1.
#' @param log.CPM.transform a logical value. If TRUE, the source data will be transformed into log-(CPM+const) before estimating the probability distributions
#' @param lib.size.params NULL or a named numerical vector containing parameters for simulating library sizes from log-normal distribution. If lib.size.params =NULL (default), then the package will fit a log-normal distribution for the library sizes in the source data to simulate new library sizes. If the user would like to specify the parameters of the log-normal distribution for the desired library sizes, then the log-mean and log-SD params of rlnorm() functions can be passed using this argument. Example, lib.size.params = c(meanlog=10, sdlog=0.2). See also ?rlnorm. 
#' @param seed an integer between 1 and 1e10 for reproducible simulation. It will be used with #set.seed() function. See also ?set.seed
#' @param verbose a logical value, if TRUE a message about the status of the simulation will be printed on the console
#' @param  ... further arguments passed to or from other methods.
#' 
#'
#'
#' @return a list of SingleCellExperiment/list objects each containing simulated counts (not normalized), smple/cell level information in colData, and gene/feature level information in rowData.
#'
#' 
#' @details This function uses a specially designed exponential family for density estimation 
#' to constructs the distribution of gene expression levels from a given real gene expression data
#' (e.g. single-cell or bulk sequencing data), and subsequently, simulates a new  
#' from the estimated distributions.  
#' 
#' For simulation of single-cell RNA-seq data (or any zero inflated gene expression data), the programm
#' involves an additional step to explicitly account for the high abundance of zero counts (if required). 
#' This step models the probability of zero counts as a function the mean 
#' expression of the gene and the library size of the cell (both in log scale) to add excess zeros. 
#' This can be done by using \emph{model.zero.prob=TRUE}. Note that,
#' for extremly large size data, it is recomended to use a random sample of cells to
#' reduce computation time. To enable this, add the argument \emph{subset.data=TRUE} and you 
#' can specify the number of cells to be used using \emph{n.samples} argument. 
#' For example \emph{n.samples=400}.
#' 
#' Given known groups of samples/cells in the source data, DGE is simulated by independently 
#' sampling data from distributions constructed for each group seprately. In particular, this procedure is 
#' applied on a set of genes with absolute log-fold-change in the source data more than a given threshold (\emph{lfc.thrld}). 
#' Moreover, when  the source dataset involves samples/cells processed in different batches, our 
#' simulation procedure incorporates this batch effect in the simulated data, if required. 
#' 
#' Different experimental designs can be simulated using the group and batch configuration arguments to
#' simulate biologica/experimental conditions and batchs, respectively. 
#' 
#' Also, it is important to filter the source data so that genes with suffient expression will be used to 
#' estimate the probability distributions.
#' 
#' 
#' @references
#' \itemize{ 
#' \item Assefa, A. T., Vandesompele, J., & Thas, O. (2019). SPsimSeq: semi-parametric simulation of bulk and single cell RNA sequencing data. \emph{bioRxiv}, doi: https://doi.org/10.1101/677740.
#' \item Efron, B., & Tibshirani, R. (1996). Using specially designed exponential families for density estimation. \emph{The Annals of Statistics}, 24(6), 2431-2461.
#' }
#' 
#' 
#' @examples 
#' \donttest{ 
#' #----------------------------------------------------------------
#' # Example 1: simulating bulk RNA-seq
#'  
#' # load the Zhang bulk RNA-seq data (availabl with the package) 
#' data("zhang.data.sub") 
#' 
#' # filter genes with sufficient expression (important step to avoid bugs) 
#' zhang.counts <- zhang.data.sub$counts[rowSums(zhang.data.sub$counts > 0)>=5, ]  
#' MYCN.status  <- zhang.data.sub$MYCN.status 
#' 

#' # We simulate only a single data (n.sim = 1) with the following property
#' # - 2000 genes ( n.genes = 2000) 
#' # - 20 samples (tot.samples = 20) 
#' # - the samples are equally divided into 2 groups each with 90 samples 
#' #   (group.config = c(0.5, 0.5))
#' # - all samples are from a single batch (batch = NULL, batch.config = 1)
#' # - we add 10% DE genes (pDE = 0.1) 
#' # - the DE genes have a log-fold-change of at least 0.5 in 
#' #   the source data (lfc.thrld = 0.5)
#' # - we do not model the zeroes separately, they are the part of density 
#' #    estimation (model.zero.prob = FALSE)
#' 
#' # simulate data
#' set.seed(6452)
#' zhang.counts2 <- zhang.counts[sample(nrow(zhang.counts), 2000), ]
#' sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = zhang.counts2, batch = NULL,
#'                           group = MYCN.status, n.genes = 2000, batch.config = 1,
#'                           group.config = c(0.5, 0.5), tot.samples = 20, 
#'                           pDE = 0.1, lfc.thrld = 0.5, 
#'                           result.format = "list",  seed = 2581988)
#' 
#' sim.data.bulk1 <- sim.data.bulk[[1]]                              
#' head(sim.data.bulk1$counts[, seq_len(5)])  # count data
#' head(sim.data.bulk1$colData)        # sample info
#' head(sim.data.bulk1$rowData)        # gene info
#' 
#' 
#' #----------------------------------------------------------------
#' # Example 2: simulating single cell RNA-seq from a single batch (read-counts)
#' # we simulate only a single scRNA-seq data (n.sim = 1) with the following property
#' # - 2000 genes (n.genes = 2000) 
#' # - 100 cells (tot.samples = 100) 
#' # - the cells are equally divided into 2 groups each with 50 cells 
#' #   (group.config = c(0.5, 0.5))
#' # - all cells are from a single batch (batch = NULL, batch.config = 1)
#' # - we add 10% DE genes (pDE = 0.1) 
#' # - the DE genes have a log-fold-change of at least 0.5
#' #' # - we model the zeroes separately (model.zero.prob = TRUE)
#' # - the ouput will be in SingleCellExperiment class object (result.format = "SCE")
#' 
#' 
#' library(SingleCellExperiment)
#' 
#' # load the NGP nutlin data (availabl with the package, processed with 
#' # SMARTer/C1 protocol, and contains read-counts)
#' data("scNGP.data")
#' 
#' # filter genes with sufficient expression (important step to avoid bugs) 
#' scNGP.data2 <- scNGP.data[rowSums(counts(scNGP.data) > 0)>=5, ]  
#' treatment <- ifelse(scNGP.data2$characteristics..treatment=="nutlin",2,1)
#' 
#' set.seed(654321)
#' scNGP.data2 <- scNGP.data2[sample(nrow(scNGP.data2), 2000), ]
#' 
#' # simulate data (we simulate here only a single data, n.sim = 1)
#' sim.data.sc <- SPsimSeq(n.sim = 1, s.data = scNGP.data2, batch = NULL,
#'                         group = treatment, n.genes = 2000, batch.config = 1,
#'                         group.config = c(0.5, 0.5), tot.samples = 100, 
#'                         pDE = 0.1, lfc.thrld = 0.5, model.zero.prob = TRUE,
#'                         result.format = "SCE",  seed = 2581988)
#' 
#' sim.data.sc1 <- sim.data.sc[[1]]
#' class(sim.data.sc1)
#' head(counts(sim.data.sc1)[, seq_len(5)])
#' colData(sim.data.sc1)
#' rowData(sim.data.sc1)
#' 
#' }
#' 
#' 
#' @export     
#' @importFrom MASS mvrnorm
#' @importFrom stats pnorm runif rbinom predict approx quantile glm rlnorm lm sd var
#' @importFrom methods is
#' @importFrom SingleCellExperiment counts colData rowData SingleCellExperiment
#' @importFrom Hmisc cut2
#' @importFrom fitdistrplus fitdist
#' @importFrom stats glm pnorm coef vcov
#' @importFrom edgeR calcNormFactors
#' @importFrom SingleCellExperiment counts SingleCellExperiment colData rowData
#' @importFrom graphics hist
#' @importFrom mvtnorm rmvnorm
#' @importFrom WGCNA cor
#' @importFrom limma voom
SPsimSeq <- function(n.sim=1, s.data, batch=NULL, group=NULL, n.genes=1000, batch.config= 1,  group.config=1, 
                     lfc.thrld=0.5, t.thrld=2.5, llStat.thrld=5, tot.samples=150, pDE=0.2, model.zero.prob=FALSE, 
                     const=1, w=NULL, geneTogenesCor=TRUE,  log.CPM.transform=TRUE, lib.size.params=NULL,
                     result.format="SCE", verbose=TRUE,  seed=2581988, ...)
{
  # Quick checks for error in the inputs
  checkInputs <- checkInputValidity(s.data = s.data, group = group, batch = batch, 
                                    group.config = group.config, batch.config = batch.config)
  
  if(any(checkInputs$val==-1)){
    msgs <- checkInputs$mes[checkInputs$val==-1]
    message(paste0("Error ", seq_len(length(msgs)), "-- ", msgs)) 
    if(length(msgs) > 1){
      stop(paste("There are", length(msgs),"errors!"),call. = FALSE)
    }else{
      stop(paste("There is", length(msgs),"error!"),call. = FALSE)
    }
  }else if(any(checkInputs$val==1)){
    msgs <- checkInputs$mes[checkInputs$val==1]
    message(paste0("Note ", seq_len(length(msgs)), "-- ", msgs)) 
    if(length(msgs) > 1){
      message(paste("There are", length(msgs)," messages."),call. = FALSE)
    }else{
      message(paste("There is", length(msgs)," message."),call. = FALSE)
    } 
  }
  
  # experiment configurartion
  if(verbose) {message("Configuring design ...")}
  if(!is.null(group)){ 
    group <- as.numeric(factor(group, labels =  seq_len(length(unique(group))))) 
  }
  if(!is.null(batch)){ 
    batch <- as.numeric(factor(batch, labels =  seq_len(length(unique(batch)))))
  }
  
  null.group = ifelse(!is.null(group), which.max(table(group))[[1]], 1)
  exprmt.design <- expriment.config(batch.config = batch.config, group.config = group.config,
                                    tot.samples = tot.samples)
  n.batch <- exprmt.design$n.batch
  n.group <- exprmt.design$n.group
  config.mat <- exprmt.design$exprmt.config 
  
  # sim control 
  if(!is.null(seed)){set.seed(seed)} 
  
  # prepare source data
  if(verbose) {message("Preparing source data ...")}
  prepare.S.Data <- prepareSourceData(s.data, batch = batch, group = group, 
                                      exprmt.design=exprmt.design, const=const,
                                      lfc.thrld=lfc.thrld, t.thrld=t.thrld, llStat.thrld=llStat.thrld, 
                                      simCtr=NULL, w=w, log.CPM.transform=log.CPM.transform,
                                      lib.size.params=lib.size.params, ...)
  LL         <- prepare.S.Data$LL
  cpm.data   <- prepare.S.Data$cpm.data
  sub.batchs <- prepare.S.Data$sub.batchs
  
  if(is.null(group)) group <- rep(1, ncol(s.data))
  if(is.null(batch)) batch <- rep(1, ncol(s.data))
  
  # fit logistic regression for the probability of zeros
  if(model.zero.prob){
    if(verbose) {message("Fitting zero probability model ...")}
    if(is(s.data, "SingleCellExperiment")){
      fracZero.logit.list <- lapply(unique(sub.batchs), function(b){
        zeroProbModel(cpm.data = cpm.data[, batch==b], L=colSums(counts(s.data)[, batch==b]),
                      simCtr=NULL, const=const, ...)
      }) 
    }
    else if(is(s.data, "data.frame") | is(s.data, "matrix")){
      fracZero.logit.list <- lapply(unique(sub.batchs), function(b){
        zeroProbModel(cpm.data = cpm.data[, batch==b], L=colSums(s.data[, batch==b]), 
                      simCtr=NULL, const=const,...)
      })
    }
  }
  
 
  # candidate genes
  if(verbose) {message("Selecting genes ...")}
  null.genes0     <- prepare.S.Data$cand.genes$null.genes
  nonnull.genes0  <- prepare.S.Data$cand.genes$nonnull.genes
  
  if(pDE>0 & !is.null(group) & length(nonnull.genes0)==0){
    warning("No gene met the criterial to be a candidiate DE gene. Perhaps consider lowering the 'lfc.thrld' or the 'llStat.thrld' or the 't.thrld'. Consequently, all the simulated genes are not DE.")
  }
  
  # Obtain copulas
  if(verbose & geneTogenesCor) {message("Estimating batch specific copulas ...")}
  if(geneTogenesCor){
    copulas.batch <- genesCopula(X = cpm.data, batch = batch, n.batch=n.batch)
  }
  else{
    copulas.batch <- lapply(sort(unique(batch)), function(bb){ 
      U    <- matrix(runif(n.batch[bb]*nrow(cpm.data), 0, 1), n.batch[bb], nrow(cpm.data))   
      colnames(U) <- rownames(cpm.data)
      U
    }) 
  }
  
  
  # simulation step
  if(verbose) {message("Simulating data ...")}
  sim.data.list <- lapply(seq_len(n.sim), function(h){
    if(verbose) {message(" ...", h, " of ", n.sim)}
    
    # sample DE and null genes 
    selctGenes <- selectGenesSim(pDE = pDE, group = group, n.genes = n.genes, null.genes0 = null.genes0,
                                 nonnull.genes0 = nonnull.genes0, group.config = group.config, h=h)
    DE.ind <- selctGenes$DE.ind
    sel.genes <- selctGenes$sel.genes
    
    
    min.val <- log(const)
    # estimate batch specific parameters
    est.list <- lapply(sel.genes, function(i){
      #print(i)
      gene.parm.est(cpm.data.i = cpm.data[i, ], batch = batch, group = group, 
                    null.group = null.group, sub.batchs = sub.batchs, de.ind = DE.ind[i], 
                    model.zero.prob = model.zero.prob, min.val = min.val, w=w, ...)
    })
    
    sim.list <- lapply(seq_len(length(sel.genes)), function(i){
      #print(i)
      if(!is.null(est.list[[i]]) & 
         !any(sapply(est.list[[i]]$batch.est, is.null)))  
      {
        # print(i) 
        if(DE.ind[i]==0){
          if(length(n.batch)>1){   
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
            par.sample <- lapply(sort(unique(group)), function(g){
     
              par.sample.g <- as.matrix(do.call("rbind.fill2", 
                                                lapply(est.list[[i]]$batch.est[[g]], function(bt){ 
                                                  data.frame(t(as.matrix(c(mvrnorm(n = 1, mu= bt$parm.list$betas, 
                                                                                   Sigma = bt$parm.list$v), 
                                                                           mu.hat=bt$parm.list$mu.hat,
                                                                           sig.hat=bt$parm.list$sig.hat)))) 
                                                })))
              
              par.sample.g 
            }) 
          }
          else if(length(n.batch)==1){
            par.sample <- lapply(sort(unique(group)), function(g){
              par.sample.g <-t(as.matrix(est.list[[i]]$Mu.batch[[g]]))
    
              par.sample.g 
            }) 
          } 
        }
        
        
        # estimate carrier density (g0)
        if(DE.ind[i]==0){
          g0 <- lapply(seq_len(nrow(par.sample)), function(bb){
            b.data <- est.list[[i]]$batch.est[[bb]] 
            gg0 <- (pnorm(b.data$yy$uls, as.matrix(par.sample)[bb, "mu.hat"],
                          as.matrix(par.sample)[bb, "sig.hat"]) -
                      pnorm(b.data$yy$lls,as.matrix(par.sample)[bb, "mu.hat"], 
                            as.matrix(par.sample)[bb, "sig.hat"]))
            
            gg0[is.nan(gg0)] <- 0
            gg0
          })
        }
        else{
          g0 <- lapply(sort(unique(group)), function(g){
            par.sample.g <- par.sample[[g]] #par.sample[[paste0("grp_", g)]]
            lapply(seq_len(nrow(par.sample.g)), function(bb){  
              b.data <- est.list[[i]]$batch.est[[bb]][[g]]
              gg0 <- (pnorm(b.data$yy$uls, par.sample.g[bb, "mu.hat"][[1]], 
                            par.sample.g[bb, "sig.hat"][[1]]) -
                        pnorm(b.data$yy$lls, par.sample.g[bb, "mu.hat"][[1]], 
                              par.sample.g[bb, "sig.hat"][[1]]))
              gg0[is.nan(gg0)] <- 0
              gg0 
            })
          }) 
        }
        
        
        
        # estimate the density g1(y)
        if(DE.ind[i]==0){
          g1 <- lapply(seq_len(nrow(par.sample)), function(bb){ 
            b.data <- est.list[[i]]$batch.est[[bb]]
            gg0 <- g0[[bb]]*sum(b.data$yy$Ny) +1
            
            s <- b.data$yy$S
            s.mat <- matrix(NA, ncol = length(coef(b.data$llm)), nrow=length(s))
            for(i in seq_len(ncol(s.mat))){
              s.mat[, i] <- s^(i-1)
            }
            gg1 <- exp(s.mat %*% par.sample[bb, seq_len(ncol(s.mat))])*gg0
            gg1 <- data.frame(gy=gg1, s=s, lls=b.data$yy$lls, uls=b.data$yy$uls)
            gg1$gy[is.infinite(gg1$gy)] <- max(gg1$gy[!is.infinite(gg1$gy)]) 
            gg1$Gy <- cumsum(gg1$gy)/sum(gg1$gy)
            gg1 
          })  
        }
        else{
          g1 <- lapply(sort(unique(group)), function(g){
            par.sample.g <- par.sample[[g]] #par.sample[[paste0("grp_", g)]]
            g0.g <- g0[[g]] #g0[[paste0("grp_", g)]]
            
            lapply(seq_len(nrow(par.sample.g)), function(bb){   
              b.data <- est.list[[i]]$batch.est[[bb]][[g]]
              gg0 <- g0.g[[bb]]*sum(b.data$yy$Ny) + 1
              
              s <- b.data$yy$S
              s.mat <- matrix(NA, ncol = length(coef(b.data$llm)), nrow=length(s))
              for(k in seq_len(ncol(s.mat))){
                s.mat[, k] <- s^(k-1)
              }
              gg1 <- exp(s.mat %*% as.matrix(par.sample.g[bb, seq_len(ncol(s.mat))]))*gg0
              gg1 <- data.frame(gy=gg1, s=s, lls=b.data$yy$lls, uls=b.data$yy$uls)
              gg1$gy[is.infinite(gg1$gy)] <- max(gg1$gy[!is.infinite(gg1$gy)]) 
              gg1$Gy <- cumsum(gg1$gy)/sum(gg1$gy)
              gg1
            }) 
          }) 
        }
        
        # simulate data  
        if(DE.ind[i]==0){ 
          Y.star <- lapply(seq_len(nrow(par.sample)), function(bb){
            Y0 <- cpm.data[sel.genes[i], (batch==bb & group==null.group)] 
            #set.seed(sim.seed)
            u <-  copulas.batch[[bb]][, sel.genes[i]]# #runif(n.batch[[bb]]) #
            gg1 <- g1[[bb]]
            #set.seed(sim.seed+1)
            y.star.b <- sapply(u, function(uu){
              yy  <- gg1$s[which.min(abs(gg1$Gy-uu))]
              difs<- diff(gg1$s)
              difs<- difs[!(is.na(difs) | is.infinite(difs) | is.nan(difs))]
              eps <- abs(mean(difs, na.rm=TRUE)) #gg1$s[2]-gg1$s[1]  
              yy  <- suppressWarnings(runif(1, yy-eps/2, yy+eps/2)) 
              yy
            })
            LL.b <- as.numeric(do.call("c", LL[[bb]]))
            if(log.CPM.transform){
              y.star.b <- round(((exp(y.star.b)-const)*LL.b)/1e6)
              y.star.b[y.star.b<0] <- 0
            } 
            
            if(model.zero.prob & mean(Y0==min.val)>0.25){
              lLL_b <- log(LL.b)
              pred.pz  <- try(predict(fracZero.logit.list[[bb]], type="response",
                                      newdata=data.frame(x1=mean(Y0), x2=lLL_b)), 
                              silent = TRUE)
              if(!is(pred.pz,"try-error")){
                #set.seed(sim.seed+2)
                drop.mlt <- sapply(pred.pz, function(p){ 
                  rbinom(1, 1, p) 
                })
              }
              else{
                drop.mlt <- 0
              } 
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
            y.star.g <- lapply(seq_len(nrow(par.sample.g)), function(bb){
              Y0.g <- cpm.data[sel.genes[i], batch==bb & group==g]
              #set.seed(sim.seed+bb)
              u <- copulas.batch[[bb]][seq_len(n.group[g]), sel.genes[i]]#runif(config.mat[bb, g])#runif(n.batch[[bb]]/length(n.group))
              gg1 <- g1.g[[bb]]
              #set.seed(sim.seed+bb+1)
              y.star.b <- sapply(u, function(uu){
                yy  <- gg1$s[which.min(abs(gg1$Gy-uu))]
                difs<- diff(gg1$s)
                difs<- difs[!(is.na(difs) | is.infinite(difs) | is.nan(difs))]
                eps <- abs(mean(difs, na.rm=TRUE)) #gg1$s[2]-gg1$s[1]  
                yy  <- suppressWarnings(runif(1, yy-eps/2, yy+eps/2))
                yy
              })
              LL.b.g <- LL[[bb]][[g]]
               
              if(log.CPM.transform){
                y.star.b <- round(((exp(y.star.b)-const)*LL.b.g)/1e6)
                y.star.b[y.star.b<0] <- 0
              }  
              
              if(model.zero.prob & mean(Y0.g==min.val)>0.25){
                lLL.b.g<- log(LL.b.g)
                pred.pz  <- try(predict(fracZero.logit.list[[bb]], type="response",
                                        newdata=data.frame(x1=mean(Y0.g), x2=lLL.b.g)), 
                                silent = TRUE)
                if(!is(pred.pz,"try-error")){
                  #set.seed(sim.seed+bb+2)
                  drop.mlt <- sapply(pred.pz, function(p){
                    ##set.seed(sim.seed+bb+2)
                    rbinom(1, 1, p)
                  })
                }
                else{
                  drop.mlt <- 0
                } 
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
    sim.count[is.na(sim.count)] <- 0
    rownames(sim.count) <- paste0("Gene_", seq_len(nrow(sim.count)))
    colnames(sim.count) <- paste0("Sample_", seq_len(ncol(sim.count)))
    
    col.data <- data.frame(Batch=rep(rep(seq_len(length(n.batch)), times=length(n.group)), config.mat),
                           Group=rep(rep(seq_len(length(n.group)), each=length(n.batch)), config.mat), 
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