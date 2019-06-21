# Sampling new data from the estimated density
SPsampleData <- function(de.ind, par.sample, cpm.data.i, batch, group, null.group,
                         g1, LL, model.zero.prob, min.val, const, config.mat, 
                         n.batch, fracZero.logit.list=NULL, ...){
  if(de.ind==0){ 
    Y.star <- lapply(1:nrow(par.sample), function(bb){
      Y0 <- cpm.data.i[(batch==bb & group==null.group)] 
      #set.seed(sim.seed)
    
        
      u <- runif(n.batch[[bb]]) 
      gg1 <- g1[[bb]]
      #smooth mode;
      #m.Gy <- loess(s~log(Gy+1), data = gg1)
      rownames(gg1) <- 1:nrow(gg1)
      #set.seed(sim.seed+1)
      y.star.b <- sapply(u, function(uu){
        whc.r  <- which.min(abs(gg1$Gy-uu))[1] # incase there are multiples
        yy     <- gg1$s[whc.r] 
        
        #yy <- predict(m.Gy, log(uu+1))
        
        
        difs<- diff(gg1$s)
        difs<- difs[!(is.na(difs) | is.infinite(difs) | is.nan(difs))]
        eps <- abs(mean(difs, na.rm=TRUE)) #gg1$s[2]-gg1$s[1]
        
        # if(whc.r==1){
        #   yll <- max(yy - 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #              gg1$lls.expd[whc.r])
        #   yul <- min(min(yy + 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #                  gg1$uls.expd[whc.r]), gg1$s[whc.r+1])
        # }else if(whc.r==nrow(gg1)){
        #   yll <- max(gg1$s[whc.r-1], max(yy - 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #                                  gg1$lls.expd[whc.r]))
        #   yul <- min(yy + 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #              gg1$uls.expd[whc.r])
        # }else{
        #   yll <- max(gg1$s[whc.r-1], max(yy - 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #                                  gg1$lls.expd[whc.r]))
        #   yul <- min(min(yy + 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
        #                  gg1$uls.expd[whc.r]), gg1$s[whc.r+1])
        # }
        yy  <- suppressWarnings(runif(1, yy-eps/2, yy+eps/2))
        yy
      })
      LL.b <- as.numeric(do.call("c", LL[[bb]]))
      y.star.b <- round(((exp(y.star.b)-const)*LL.b)/1e6)
      y.star.b[y.star.b<0] <- 0
      
      if(model.zero.prob & mean(Y0==min.val)>0.25){
        lLL_b <- log(LL.b)
        pred.pz  <- try(predict(fracZero.logit.list[[bb]], type="response",
                                newdata=data.frame(x1=mean(Y0), x2=lLL_b)), 
                        silent = TRUE)
        if(class(pred.pz) != "try-error"){
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
      y.star.g <- lapply(1:nrow(par.sample.g), function(bb){
        Y0.g <- cpm.data.i[batch==bb & group==g]
        #set.seed(sim.seed+bb)
        u <- runif(config.mat[bb, g])#runif(n.batch[[bb]]/length(n.group))
        gg1 <- g1.g[[bb]]
        rownames(gg1) <- 1:nrow(gg1)
        #smooth mode;
        #m.Gy.g <- loess(s~log(Gy+1), data = gg1)
        #set.seed(sim.seed+bb+1)
        y.star.b <- sapply(u, function(uu){
          whc.r  <- which.min(abs(gg1$Gy-uu))[1] # incase there are multiples
          yy     <- gg1$s[whc.r]
          
          #yy <- predict(m.Gy.g, log(uu+1))
          difs<- diff(gg1$s)
          difs<- difs[!(is.na(difs) | is.infinite(difs) | is.nan(difs))]
          eps <- abs(mean(difs, na.rm=TRUE)) #gg1$s[2]-gg1$s[1]
          # if(whc.r==1){
          #   yll <- max(yy - 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
          #              gg1$lls.expd[whc.r])
          #   yul <- min(min(yy + 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
          #                  gg1$uls.expd[whc.r]), gg1$s[whc.r+1])
          # }else if(whc.r==nrow(gg1)){
          #   yll <- max(gg1$s[whc.r-1], max(yy - 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
          #                                  gg1$lls.expd[whc.r]))
          #   yul <- min(yy + 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
          #              gg1$uls.expd[whc.r])
          # }else{
          #   yll <- max(gg1$s[whc.r-1], max(yy - 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
          #                                  gg1$lls.expd[whc.r]))
          #   yul <- min(min(yy + 0.25*abs(gg1$uls.expd[whc.r]-gg1$lls.expd[whc.r]),
          #                  gg1$uls.expd[whc.r]), gg1$s[whc.r+1])
          # }
          yy  <- suppressWarnings(runif(1, yy-eps/2, yy+eps/2))
          
          yy
        })
        LL.b.g <- LL[[bb]][[g]]
        y.star.b <- round(((exp(y.star.b)-const)*LL.b.g)/1e6)
        y.star.b[y.star.b<0] <- 0
        
        if(model.zero.prob & mean(Y0.g==min.val)>0.25){
          lLL.b.g<- log(LL.b.g)
          pred.pz  <- try(predict(fracZero.logit.list[[bb]], type="response",
                                  newdata=data.frame(x1=mean(Y0.g), x2=lLL.b.g)), 
                          silent = TRUE)
          if(class(pred.pz) != "try-error"){
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
  }
  Y.star
}






