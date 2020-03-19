#'Select candidate genes
#'
#'This function can be used to independently select candidate genes from a given real RNA-srq data (bulk/single)
#'for the SPsimSeq simulation. It chooses genes with various chracteristics, such as log-fold-change 
#' above a certain thereshold. 
#' 
#' @param cpm.data logCPM transformed matrix (if log.CPM.transform=FALSE, then it is the source gene expression data)
#' @param X a vector of indicators for group memebership of cells/samples 
#' @param lfc.thrld a positive numeric value for the minimum absolute log-fold-change for selecting candidate DE genes in the source data (when group is not NULL and pDE>0)
#' @param t.thrld a positive numeric value for the minimum absolute t-test statistic for the log-fold-changes of genes for selecting candidate DE genes in the source data (when group is not NULL and  pDE>0)
#' @param llStat.thrld a positive numeric value for the minimum squared test statistics from the log-linear model to select candidate DE genes in the source data (when group is not NULL and  pDE>0)
#' containing X as a covariate to select DE genes 
#' @param max.frac.zeror.diff a numeric value >=0 indicating the maximum absolute
#' difference in the fraction of zero counts between the groups for DE genes.
#' @param w a numeric value between 0 and 1. The number of classes to construct the probability distribution will be round(w*n), where n is the total number of samples/cells in a particular batch of the source data
#' @param const a positive constant to be added to the CPM before log transformation, to avoid log(0). The default is 1.
#' @return a list object contating a set of candidate null and non-null genes and additional results
#' @examples
#'  # example: see ?SPsimSeq
  
#' @importFrom stats lm sd density rnbinom rlnorm var runif predict rbinom rgamma
#' @importFrom utils combn
chooseCandGenes <- function(cpm.data, X,  lfc.thrld, llStat.thrld, t.thrld,
                             max.frac.zeror.diff = Inf, const, w){
  n.cells   <- table(X)
  sim.group <- length(n.cells)
  # calculate fold-changes
  logConst =  log(const)
  m.diff  <- apply(cpm.data, 1, function(y){ 
    l.mod  <- lm(y~X)
    t.stat <- max(abs(summary(l.mod)[["coefficients"]][-1, "t value"]))
    #Removed as.numeric
    fc     <- max(abs(coef(l.mod)[-1]))
    frac.z.diff <- max(abs(combn(tapply(y, X, function(yy) mean(yy==logConst)), 2, FUN=diff)))
    c(t.stat = t.stat, fc = fc, frac.z.diff = frac.z.diff)
  })
  
  null.genes0    <- colnames(m.diff)[m.diff["t.stat",] < t.thrld | m.diff["fc",] < lfc.thrld ]
  nonnull.genes0 <- colnames(m.diff)[(m.diff["t.stat",] >= t.thrld) & 
                                       (m.diff["fc",] >= lfc.thrld) & 
                                       (m.diff["frac.z.diff",] <= max.frac.zeror.diff)]

  if(llStat.thrld > 0 & length(nonnull.genes0)>=1){
    statLLmodel <- sapply(nonnull.genes0, function(j){
      Y <- lapply(names(n.cells), function(x){
        cpm.data[j, X==x]
      })
      
      S.list <- lapply(Y, FUN = obtCount, w=w)
      ss     <- lapply(S.list, function(x) x$S)
      lls    <- lapply(S.list, function(x) x$lls)
      uls    <- lapply(S.list, function(x) x$uls)
      Ny     <- lapply(S.list, function(x) x$Ny)
      muHats     <- vapply(Y, mean, FUN.VALUE = 0)
      sigHats     <- vapply(Y, sd, FUN.VALUE = 0)
      N = vapply(Ny, sum, FUN.VALUE = 0)
      
      gg0 <- lapply(seq_along(Y), function(l){ 
        (pnorm(uls[[l]], muHats[[l]], sigHats[[l]]) - 
           pnorm(lls[[l]], muHats[[l]], sigHats[[l]]))*N[[l]]
      })
      Xy <- lapply(seq_along(Y), function(l) rep(l-1, length(ss[[l]])))
      gg0 = do.call("c", gg0)
      Ny <- do.call('c', Ny)
      ss <- do.call('c', ss)
      Xy <- do.call('c', Xy)
      ofs = log(gg0+1)
      
      formulae = paste0("Ny~", 
                        c("I(ss)+ I(ss^2)+ I(Xy) + I(ss*Xy) + I((ss^2)*Xy", 
                   "I(ss)+ I(ss^2)+ I(Xy) + I(ss*Xy)",
                   "I(ss)+ I(Xy) + I(ss*Xy)"
                   ))
      for (form in formulae){
      l.mod.x <- tryCatch(glm(form), 
                              family = "poisson", offset = ofs,
                          error=function(e){}, warning=function(w){})
      if(!is.null(l.mod.x)) break
      }
      if(!is.null(l.mod.x)){
        if(l.mod.x$rank == ncol(l.mod.x$R)){
          Z.X <- summary(l.mod.x)$coefficients[, 3]
          Z.X <- Z.X[names(Z.X) %in% c("I(Xy)", "I(ss * Xy)", "I((ss^2) * Xy)")]
          sum.square.Z.X <- sum(Z.X^2, na.rm = TRUE)
          sum.square.Z.X
        }
        else{0} 
      }
      else {0} 
    })
  } else if(llStat.thrld==0 & length(nonnull.genes0)>=1){
    statLLmodel <- rep(0, length(nonnull.genes0))
  } else{
    statLLmodel <- 0
    warning("Unable to find candidate non-null genes with. Consider to lower the fold-change threshold or llStat threshold.")
  }
  compr.stat2  <- t(m.diff[,nonnull.genes0])
  compr.stat2$statLL <- statLLmodel[rownames(compr.stat2)]
  nonnull.genes <- nonnull.genes0[statLLmodel>=llStat.thrld]
  null.genes    <- c(null.genes0, setdiff(nonnull.genes0, nonnull.genes))
  sel.genes <- list(null.genes=unique(null.genes) , nonnull.genes=nonnull.genes,
                    statLLmodel=statLLmodel, compr.stat=compr.stat2)
  sel.genes
}
