#' Select candidate genes
#'
#' This function can be used independently to select candidate genes from a given real RNA-srq data (bulk/single)
#' for the SPsimSeq simulation. It chooses genes with cetrail chracteristics, such as log-fold-change 
#' above a certain thereshold. 
#' 
#' @param cpm.data CPM matrix
#' @param X a vector of indicators for group memebership of cells/samples 
#' @param lfc.thrld a numeric value for the minimum fold change for DE genes
#' @param t.thrld a numeric value for the minimum t statistic for DE genes 
#' @param llStat.thrld a numeric value for the minimum squared test statistics from a log-linear model
#' containing X as a covariate to select DE genes
#' @param carrier.dist a character indicating the type of
#' carrier density (carrier.dist="normal" or carrier.dist="kernel") 
#' @param max.frac.zeror.diff a numeric value >=0 indicating the maximum  absolute
#' difference in the fraction of zero counts between the groups for DE genes. Default in Inf
#' @param const a small constant (>0) added to the CPM before log transformation, to avoid  log(0).
#' default is 1e-5, see in 
#' @param  ... further arguments passed to or from other methods.
#'
#' @return a list object contating a set of candidate null and non-null genes and additional results
#' @examples
#'  # example
#' @export  
#' @importFrom stats lm sd density rnbinom rlnorm var runif predict rbinom rgamma
#' @importFrom SingleCellExperiment counts
chooseCandGenes <- function(cpm.data, X,  lfc.thrld=1,  
                             llStat.thrld=10, t.thrld=2.5, carrier.dist="normal",
                             max.frac.zeror.diff=Inf, const=1e-5,...){
  n.cells   <- table(X)
  sim.group <- length(n.cells)
  
  # calculate log CPM
  # if(class(s.data)=="SingleCellExperiment"){
  #   cpm.data <- log(calCPM(counts(s.data))+1) 
  # }
  # else if(class(s.data) %in% c("data.frame", "matrix")){
  #   cpm.data <- log(calCPM(s.data)+1) 
  # } 

  # calculate fold-changes
  m.diff  <- as.data.frame(t(apply(cpm.data, 1, function(y){ 
    l.mod  <- lm(y~X)
    t.stat <- max(abs(as.numeric(summary(l.mod)[["coefficients"]][-1, "t value"])))
    fc     <- max(abs(as.numeric(coef(l.mod)[-1])))
    frac.z.diff <- max(abs(combn2(tapply(y, X, function(yy) mean(yy==log(const))), 2, FUN=diff)))
    c(t.stat=t.stat, fc=fc, frac.z.diff=frac.z.diff)
  })))
  
  null.genes0      <- rownames(m.diff)[m.diff$t.stat<t.thrld | m.diff$fc<lfc.thrld ]
  nonnull.genes0   <- rownames(m.diff)[(m.diff$t.stat>=t.thrld) & (m.diff$fc>=lfc.thrld) & (m.diff$frac.z.diff <= max.frac.zeror.diff)]
  
  compr.stat <- m.diff 
  statLLmodel <- sapply(nonnull.genes0, function(j){
    Y <- lapply(names(n.cells), function(x){
      as.numeric(cpm.data[j, X==x])
    })
    
    S.list <- lapply(Y, obtCount)
    ss     <- lapply(S.list, function(x) x$S)
    lls    <- lapply(S.list, function(x) x$lls)
    uls    <- lapply(S.list, function(x) x$uls)
    Ny     <- lapply(S.list, function(x) x$Ny)
    w      <- sapply(S.list, function(x) x$w)[[1]]

    N=sapply(Ny, sum)

    if(carrier.dist=="kernel"){
      g0 <- lapply(1:length(Y), function(l){
        density(Y[[l]], from=min(ss[[l]])-w/2, to=max(ss[[l]])+w/2)
      })
      gg0 <- lapply(1:length(Y), function(l){
        sapply(ss[[l]], function(s) g0[[l]]$y[which.min(abs(g0[[l]]$x-s))])*N[[l]]
      })
    }
    else if(carrier.dist=="normal"){
      gg0 <- lapply(1:length(Y), function(l){ 
        mu.hat <- mean(Y[[l]])
        sig.hat<- sd(Y[[l]]) 
        (pnorm(uls[[l]], mu.hat, sig.hat)-pnorm(lls[[l]], mu.hat, sig.hat))*N[[l]]
      })
    }
     
    Xy <- lapply(1:length(Y), function(l) rep(l-1, length(ss[[l]])))

    ofs <- 1

    Ny <- do.call('c', Ny)
    gg0<- do.call('c', gg0)
    ss <- do.call('c', ss)
    Xy <- do.call('c', Xy)


    # l.mod.x <- try(glm(Ny~I(ss)+ I(ss^2)+ I(Xy) + I(ss*Xy) + I((ss^2)*Xy),
    #                    family = "poisson", offset = log(gg0+ofs)),
    #                silent = TRUE)
    # if(all(class(l.mod.x) != "try-error")){
    #   if(l.mod.x$rank != ncol(l.mod.x$R)){
    #     l.mod.x <- try(glm(Ny~I(ss)+ I(ss^2)+ I(Xy) + I(ss*Xy),
    #                        family = "poisson", offset = log(gg0+ofs)),
    #                    silent = TRUE)
    #     if(all(class(l.mod.x) != "try-error")){
    #       if(l.mod.x$rank != ncol(l.mod.x$R)){
    #         l.mod.x <- try(glm(Ny~I(ss)+ I(Xy) + I(ss*Xy),
    #                            family = "poisson", offset = log(gg0+ofs)),
    #                        silent = TRUE)
    #       } 
    #     }
    #   } 
    # }
    
    l.mod.x <- tryCatch(glm(Ny~I(ss)+ I(ss^2)+ I(Xy) + I(ss*Xy) + I((ss^2)*Xy), 
                                   family = "poisson", offset = log(gg0+ofs)),
                    error=function(e){}, warning=function(w){})
    if(is.null(l.mod.x)){
      l.mod.x <- tryCatch(glm(Ny~I(ss)+ I(ss^2)+ I(Xy) + I(ss*Xy),
                                     family = "poisson", offset = log(gg0+ofs)),
                      error=function(e){}, warning=function(w){})
      if(is.null(l.mod.x)){
        l.mod.x <- tryCatch(suppressWarnings(glm(Ny~I(ss)+ I(Xy) + I(ss*Xy),
                                                       family = "poisson", 
                                                       offset = log(gg0+ofs))),
                        error=function(e){})
        if(is.null(l.mod.x)){
          l.mod.x <- NULL
        } 
      } 
    } 
 
    if(!is.null(l.mod.x)){
      #coef.X <- coef(l.mod.x)[grep("Xy", names(coef(l.mod.x)))]
      #sum.square.coef.X <- sum(coef.X^2, na.rm = TRUE)
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

  compr.stat2  <- compr.stat[nonnull.genes0, ]
  compr.stat2$statLL <- statLLmodel[rownames(compr.stat2)]
  top.genes0 <- nonnull.genes0[statLLmodel>=llStat.thrld]

  nonnull.genes <- top.genes0
  null.genes    <- c(null.genes0, nonnull.genes0[!(nonnull.genes0 %in% top.genes0)])
  sel.genes <- list(null.genes=unique(null.genes) , nonnull.genes=nonnull.genes,
                    statLLmodel=statLLmodel, compr.stat=compr.stat2)
  sel.genes

}
