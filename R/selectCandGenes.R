#' select candidate genes
#'
#' More detailed description
#'
#' @param s.data a source data (a SingleCellExperiment object)
#' @param X a vector of indicators for group memebership of cells/samples 
#' @param lfc.thrld a numeric value for the minimum fold change for DE genes
#' @param t.thrld a numeric value for the minimum t statistic for DE genes
#' ('mean.diff'= difference in the mean log CPM, 'rank.sum'= a U-statistic for the log CPM)
#' @param llStat.thrld a numeric value for the minimum squared test statistics from a log-linear model
#' containing X as a covariate to select DE genes
#' @param carrier.dist a character indicating the type of
#' carrier density (carrier.dist="normal" or carrier.dist="kernel")
#' @param max.frac.zero a numeric value between 0 and 1 indicating the maximum fraction of
#'  zero counts that a DE gene should have
#' @param max.frac.zeror.diff a numeric value between 0 and 1 indicating the maximum  absolute
#' difference in the fraction of zero counts between the groups for DE genes
#' @param  ... further arguments passed to or from other methods.
#'
#' @return a list object contating a set of candidate DE and null genes and additional results
#' @examples
#'  # example
#' @export 
#' @importFrom utils combn

chooseCandGenes <- function(s.data, X,  lfc.thrld=0,  
                             llStat.thrld=10, t.thrld=2.5, carrier.dist="normal",
                             max.frac.zero=1, max.frac.zeror.diff=Inf, ...){
  n.cells   <- table(X)
  sim.group <- length(n.cells)
  
  # calculate log CPM
  if(class(s.data)=="SingleCellExperiment"){
    cpm.data <- log(calCPM(counts(s.data))+1) 
  }
  else if(class(s.data) %in% c("data.frame", "matrix")){
    cpm.data <- log(calCPM(s.data)+1) 
  } 

  # calculate fold-changes
  m.diff  <- as.data.frame(t(apply(cpm.data, 1, function(y){ 
    l.mod  <- lm(y~X)
    t.stat <- max(abs(as.numeric(summary(l.mod)[["coefficients"]][-1, "t value"])))
    fc     <- max(abs(as.numeric(coef(l.mod)[-1])))
    frac.z.diff <- max(abs(combn(tapply(y, X, function(yy) mean(yy==0)), 2, FUN=diff)))
    c(t.stat=t.stat, fc=fc, frac.z.diff=frac.z.diff)
  })))
  
  null.genes0      <- rownames(m.diff)[m.diff$t.stat<t.thrld | m.diff$fc<lfc.thrld ]
  nonnull.genes0   <- rownames(m.diff)[(m.diff$t.stat>=t.thrld) & (m.diff$fc>=lfc.thrld) & (m.diff$frac.z.diff <= max.frac.zeror.diff)]
  
  compr.stat <- m.diff 
  statLLmodel <- sapply(nonnull.genes0, function(j){
    Y <- lapply(names(n.cells), function(x){
      as.numeric(cpm.data[j, X==x])
    })

    # s <- lapply(1:length(Y), function(l){
    #   if(!is.null(w)){
    #     ww=w
    #     while(round(ww*length(Y[[l]]))<3 & ww<1){
    #       ww <- ww+0.05
    #     }
    #     h <- hist(Y[[l]], nclass = round(ww*length(Y[[l]])), plot = FALSE, right = TRUE)
    #   }
    #   else{
    #     h <- hist(Y[[l]], plot = FALSE, right = TRUE)
    #   }
    # 
    #   h$breaks
    # })
    # 
    # lls <- lapply(1:length(Y), function(l){
    #   t=s[[l]]
    #   t[1:(length(t)-1)]
    # })
    # uls <- lapply(1:length(Y), function(l){
    #   t=s[[l]]
    #   t[2:length(t)]
    # })
    # ss <- lapply(1:length(Y), function(l){
    #   t1 <- lls[[l]]
    #   t2 <- uls[[l]]
    #   (t1+t2)/2
    # })
    # Ny <- lapply(1:length(Y), function(l){
    #   t  <- ss[[l]]
    #   t1 <- lls[[l]]
    #   t2 <- uls[[l]]
    #   sapply(1:length(t), function(x){
    #     if(x==1) sum(Y[[l]]<=t2[x])
    #     else sum(Y[[l]]>t1[x] & Y[[l]]<=t2[x])
    #   })
    # }) 
    
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
        # est.parms <- try(fitdist(Y[[l]], distr = "norm", method = "mle"), silent = TRUE)
        # if(class(est.parms)=="try-error"){
        #   est.parms <- fitdist(Y[[l]], distr = "norm", method = "mme")
        # }
        # est.parms <- est.parms$estimate
        # mu.hat <- est.parms[["mean"]]
        # sig.hat<- est.parms[["sd"]]
        mu.hat <- mean(Y[[l]])
        sig.hat<- sd(Y[[l]]) 
        (pnorm(uls[[l]], mu.hat, sig.hat)-pnorm(lls[[l]], mu.hat, sig.hat))*N[[l]]
      })
    }

    # plot(ss[[1]], Ny[[1]], type="l") ; lines(ss[[2]], Ny[[2]], type="l", col=2)
    # lines(ss[[1]], gg0[[1]], type="l", lty=3, lwd=3) ; lines(ss[[2]], gg0[[2]], type="l", col=2, lty=3, lwd=3)
    #

    Xy <- lapply(1:length(Y), function(l) rep(l-1, length(ss[[l]])))

    # ss  <- lapply(1:length(Y),  function(l) ss[[l]][Ny[[l]]>0] )
    # gg0 <- lapply(1:length(Y), function(l) gg0[[l]][Ny[[l]]>0] )
    # Xy  <- lapply(1:length(Y),  function(l) Xy[[l]][Ny[[l]]>0] )
    # Ny  <- lapply(1:length(Y),  function(l) Ny[[l]][Ny[[l]]>0] )


    ofs <- 1

    Ny <- do.call('c', Ny)
    gg0<- do.call('c', gg0)
    ss <- do.call('c', ss)
    Xy <- do.call('c', Xy)


    l.mod.x <- try(glm(Ny~I(ss)+ I(ss^2)+ I(Xy) + I(ss*Xy) + I((ss^2)*Xy),
                       family = "poisson", offset = log(gg0+ofs)),
                   silent = TRUE)
    if(all(class(l.mod.x) != "try-error")){
      if(l.mod.x$rank != ncol(l.mod.x$R)){
        l.mod.x <- try(glm(Ny~I(ss)+ I(ss^2)+ I(Xy) + I(ss*Xy),
                           family = "poisson", offset = log(gg0+ofs)),
                       silent = TRUE)
        if(all(class(l.mod.x) != "try-error")){
          if(l.mod.x$rank != ncol(l.mod.x$R)){
            l.mod.x <- try(glm(Ny~I(ss)+ I(Xy) + I(ss*Xy),
                               family = "poisson", offset = log(gg0+ofs)),
                           silent = TRUE)
          } 
        }
      } 
    }

    # pred1 <- predict(l.mod.x, type="response", newdata = data.frame(ss=ss[Xy==0], Xy=0, gg0=gg0[Xy==0]))
    # pred2 <- predict(l.mod.x, type="response", newdata = data.frame(ss=ss[Xy==1], Xy=1, gg0=gg0[Xy==1]))
    #
    # lines(ss[Xy==0], pred1, type="l", lty=1, lwd=3) ; lines(ss[Xy==1], pred2, type="l", col=2, lty=1, lwd=3)
    #
    if(all(class(l.mod.x) != "try-error" )){
      #coef.X <- coef(l.mod.x)[grep("Xy", names(coef(l.mod.x)))]
      #sum.square.coef.X <- sum(coef.X^2, na.rm = TRUE)
      if(l.mod.x$rank == ncol(l.mod.x$R)){
        Z.X <- summary(l.mod.x)$coefficients[, 3]
        Z.X <- Z.X[names(Z.X) %in% c("I(Xy)", "I(ss * Xy)", "I((ss^2) * Xy)")]
        sum.square.Z.X <- sum(Z.X^2, na.rm = TRUE)
        sum.square.Z.X
      }
      else{
        0
      }

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
