#' select candidate genes
#'
#' More detailed description
#'
#' @param s.data a source data (a SingleCellExperiment object)
#' @param X a vector of indicators for group memebership of cells/samples
#' @param fc.type a character indicating how DE genes should be calculated
#' ('mean.diff'= difference in the mean log CPM, 'rank.sum'= a U-statistic for the log CPM)
#' @param lfc.thrld a numeric value for the minimum fold change for DE genes (if fc.type='mean.diff')
#' @param t.thrld a numeric value for the minimum t statistic for DE genes (if fc.type='mean.diff')
#' @param U.thrld a numeric value for the minimum U-statistic for DE genes (if fc.type='rank.sum')
#' @param llStat.thrld a numeric value for the minimum squared test statistics from a log-linear model
#' containing X as a covariate to select DE genes
#' @param carrier.dist a character indicating the type of
#' carrier density (carrier.dist="normal" or carrier.dist="kernel")
#' @param w a numeric value between 0 and 1 or NULL refering the number of classes to be created
#' for the outcome data (if NULL the algorithm to calculate breakes in graphics::hist() function will be used)
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


chooseCandGenes <- function(s.data, X, fc.type="mean.diff", lfc.thrld=0, U.thrld=0.7,
                             llStat.thrld=10, t.thrld=2.5, carrier.dist="normal",
                             w=0.5, max.frac.zero=0.7, max.frac.zeror.diff=0.1, ...){
  n.cells <- table(X)
  sim.group=length(n.cells)
  cpm.data <- log(edgeR::cpm(counts(s.data))+1)
  #cpm.data[counts(s.data)==0] <- 0

  # calculate fold-changes
  if(fc.type=="rank.sum"){
    U.stat  <- apply(cpm.data, 1, function(y){
      wlx.test <- wilcox.test(y~X)
      #wlx.test$statistic*(wlx.test$p.value <0.05)
      U <- wlx.test$statistic/prod(n.cells)
      U <- ifelse(U<0.5, 1-U, U)
      U
    })

    null.genes0      <- names(U.stat)[U.stat<U.thrld]
    nonnull.genes0   <- names(U.stat)[U.stat>=U.thrld]

    compr.stat <- as.data.frame(U.stat)
  }
  else if(fc.type=="mean.diff"){
    m.diff  <- as.data.frame(t(apply(cpm.data, 1, function(y){
      t.stat <- as.numeric(abs(t.test(y~X)$statistic))
      fc     <- as.numeric(abs(diff(tapply(y, X, mean))))
      frac.z.diff <- abs(diff(tapply(y, X, function(yy) mean(yy==0))))
      c(t.stat=t.stat, fc=fc, frac.z.diff=frac.z.diff)
    })))

    null.genes0      <- rownames(m.diff)[m.diff$t.stat<t.thrld | m.diff$fc<lfc.thrld ]
    nonnull.genes0   <- rownames(m.diff)[(m.diff$t.stat>=t.thrld) & (m.diff$fc>=lfc.thrld) & (m.diff$frac.z.diff <= max.frac.zeror.diff)]

    compr.stat <- m.diff
  }

  statLLmodel <- sapply(nonnull.genes0, function(j){
    Y <- lapply(names(n.cells), function(x){
      as.numeric(cpm.data[j, X==x])
    })

    s <- lapply(1:length(Y), function(l){
      if(!is.null(w)){
        ww=w
        while(round(ww*length(Y[[l]]))<3 & ww<1){
          ww <- ww+0.05
        }
        h <- hist(Y[[l]], nclass = round(ww*length(Y[[l]])), plot = FALSE, right = TRUE)
      }
      else{
        h <- hist(Y[[l]], plot = FALSE, right = TRUE)
      }

      h$breaks
    })

    lls <- lapply(1:length(Y), function(l){
      t=s[[l]]
      t[1:(length(t)-1)]
    })
    uls <- lapply(1:length(Y), function(l){
      t=s[[l]]
      t[2:length(t)]
    })
    ss <- lapply(1:length(Y), function(l){
      t1 <- lls[[l]]
      t2 <- uls[[l]]
      (t1+t2)/2
    })
    Ny <- lapply(1:length(Y), function(l){
      t  <- ss[[l]]
      t1 <- lls[[l]]
      t2 <- uls[[l]]
      sapply(1:length(t), function(x){
        if(x==1) sum(Y[[l]]<=t2[x])
        else sum(Y[[l]]>t1[x] & Y[[l]]<=t2[x])
      })
    })

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
        est.parms <- try(fitdist(Y[[l]], distr = "norm", method = "mle"), silent = TRUE)
        if(class(est.parms)=="try-error"){
          est.parms <- fitdist(Y[[l]], distr = "norm", method = "mme")
        }
        est.parms <- est.parms$estimate
        mu.hat <- est.parms[["mean"]]
        sig.hat<- est.parms[["sd"]]
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
    if(all(class(l.mod.x) != "try-error" ) & l.mod.x$rank != ncol(l.mod.x$R)){
      l.mod.x <- try(glm(Ny~I(ss)+ I(ss^2)+ I(Xy) + I(ss*Xy),
                         family = "poisson", offset = log(gg0+ofs)),
                     silent = TRUE)
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
