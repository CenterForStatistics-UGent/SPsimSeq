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
#' @param prior.count a positive constant to be added to the CPM before log transformation, to avoid log(0). The default is 1.
#' @return a list object contating a set of candidate null and non-null genes and additional results
#' @importFrom stats lm sd density rnbinom rlnorm var runif predict rbinom rgamma
#' @importFrom utils combn
chooseCandGenes <- function(cpm.data, group, lfc.thrld, llStat.thrld, t.thrld, w =w,
                             max.frac.zeror.diff = Inf, pDE, n.genes, prior.count){
  # calculate fold-changes and t-statistics
  logConst =  log(prior.count)
  m.diff  <- apply(cpm.data, 1, function(y){ 
    l.mod  <- lm(y~group)
    t.stat <- max(abs(summary(l.mod)[["coefficients"]][-1, "t value"]))
    #Removed as.numeric
    fc     <- max(abs(coef(l.mod)[-1]))
    frac.z.diff <- max(abs(combn(tapply(y, group, function(yy) mean(yy==logConst)), 2, FUN=diff)))
    c(t.stat = t.stat, fc = fc, frac.z.diff = frac.z.diff)
  })
  null.genes0    <- colnames(m.diff)[m.diff["t.stat",] < t.thrld | m.diff["fc",] < lfc.thrld ]
  nonnull.genes0 <- colnames(m.diff)[(m.diff["t.stat",] >= t.thrld) & 
                                       (m.diff["fc",] >= lfc.thrld) & 
                                       (m.diff["frac.z.diff",] <= max.frac.zeror.diff)]
  #For the ones exceeding the threshold,
  statLLmodel <- if(length(nonnull.genes0)>=1){
     vapply(nonnull.genes0, FUN.VALUE = numeric(1), function(j){
      fits = tapply(cpm.data[j,], group, function(Y){
        mu.hat = mean(Y)
        sig.hat = sd(Y)
        #Bin counts
        countY <- obtCount(Y = Y, w = w)
        #Fit exponential density
        llModel = fitLLmodel(yy = countY, mu.hat = mu.hat, sig.hat = sig.hat, 
                             n = length(Y))
        #define offset of carrier density
        ofs = log(llModel$g0*length(Y))
        list(ofs = ofs, llModel = llModel)
      })
      #Extract counts, offsets and midpoints, and concatenate
      df = data.frame(
        counts = unlist(lapply(fits, function(x) x$llModel$counts)),
        mids = unlist(lapply(fits, function(x) x$llModel$mids)),
        ofs = unlist(lapply(fits, function(x) x$ofs)),
        groups = factor(unlist(lapply(seq_along(fits), function(l) {
          rep(l, length(fits[[l]]$ofs))})))
      )
      l.mod.x = NULL; i = 5
      #The different models to try
      terms = c("mids", "groups", "groups:mids", "I(mids^2)", "groups:I(mids^2)")
      while (is.null(l.mod.x) && i>=3){
        form = paste("counts ~", paste(terms[seq_len(i)], collapse ="+"), "+offset(ofs)")
        l.mod.x <- tryCatch(glm(form, family = "poisson", data = df),
                            error=function(e){}, warning=function(w){})
        i = i-1
      }
      #Return squared sum of test statistics
      sum.square.Z.X = if(!is.null(l.mod.x) && l.mod.x$rank == ncol(l.mod.x$R)){
          Z.X <- summary(l.mod.x)$coefficients[, 3]
          Z.X <- Z.X[grepl(names(Z.X), pattern = "groups")]
          sum(Z.X^2, na.rm = TRUE)
      } else{0} 
      return(sum.square.Z.X)
    })
  } else{
    warning("Unable to find candidate non-null genes. Consider to lower the fold-change threshold or llStat threshold.")
    0
  }
  nonnull.genes <- nonnull.genes0[statLLmodel>=llStat.thrld]
  null.genes    <- c(null.genes0, setdiff(nonnull.genes0, nonnull.genes))
  #Throw warnings where needed
  if((1-pDE)*n.genes > length(null.genes)){
    message("Note: The number of null genes (not DE) in the source data is ",
            length(null.genes),
            " and the number of null genes required to be included in the simulated data is ", 
            round((1-pDE)*(n.genes)),
            ". Therefore, candidiate null genes are sampled with replacement.")
  }
  if(pDE*n.genes > length(nonnull.genes)){
    message("Note: The number of DE genes detected in the source data is ",
            length(nonnull.genes),
            " and the number of DE genes required to be included in the simulated data is ",
            round(pDE*n.genes),
            ". Therefore, candidiate DE genes are sampled with replacement.")
  }
  if(pDE>0 & !is.null(group) & length(nonnull.genes)==0){
    warning("No gene met the criterion to be a candidiate DE gene. Perhaps consider
            lowering the 'lfc.thrld' or the 'llStat.thrld' or the 't.thrld'. Consequently,
            all the simulated genes are not DE.")
  }
  list(null.genes = unique(null.genes), nonnull.genes = nonnull.genes)
}
