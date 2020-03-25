# A function to simulate data from the estimated density (gene level)
sampleDatSPDens <- function(cpm.data, sel.genes.i, par.sample, DE.ind.ii, null.group, config.mat,
                            copulas.batch, LL, group, batch, n.group, n.batch, g1, log.CPM.transform, 
                            const, model.zero.prob, min.val, fracZero.logit.list){
  if(DE.ind.ii==0){ 
    Y.star <- lapply(seq_len(nrow(par.sample)), function(bb){
      Y0 <- cpm.data[sel.genes.i, (batch==bb & group==null.group)] 
      #set.seed(sim.seed)
      u <-  copulas.batch[[bb]][, sel.genes.i]# #runif(n.batch[[bb]]) #
      y.star.b = makeYstarB(u = u, gg1 = g1[[bb]])
      #set.seed(sim.seed+1)
      if(log.CPM.transform){
        LL.b <- as.numeric(do.call("c", LL[[bb]]))
        y.star.b <- round(((exp(y.star.b)-const)*LL.b)/1e6)
        y.star.b[y.star.b<0] <- 0
      }else{
        LL.b <- 1
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
        }else{
          drop.mlt <- 0
        } 
        y.star.b <- y.star.b*(1-drop.mlt)
      }
      data.frame("y.new"=as.numeric(y.star.b), "Batch"=rep(names(n.batch)[bb], n.batch[bb]),
                 "Group"=rep(names(n.group), config.mat[bb, ]))
    })
    Y.star <- do.call("rbind", Y.star)
    #if(any(is.na(Y.star))){print(i)}
  }
  else{ 
    Y.star <- lapply(sort(unique(group)), function(g){
      par.sample.g <- par.sample[[g]] #par.sample[[paste0("grp_", g)]]
      g1.g <- g1[[g]] #g1[[paste0("grp_", g)]]
      y.star.g <- lapply(seq_len(nrow(par.sample.g)), function(bb){
        Y0.g <- cpm.data[sel.genes.i, batch==bb & group==g]
        #set.seed(sim.seed+bb)
        u <- copulas.batch[[bb]][seq_len(config.mat[bb, g]), sel.genes.i]#runif(config.mat[bb, g])#runif(n.batch[[bb]]/length(n.group))
        y.star.b = makeYstarB(u, g1.g[[bb]]) 
        if(log.CPM.transform){
          LL.b.g <- LL[[bb]][[g]]
          y.star.b <- round(((exp(y.star.b)-const)*LL.b.g)/1e6)
          y.star.b[y.star.b<0] <- 0
        }else{
          LL.b.g <- 1
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
          }else{
            drop.mlt <- 0
          } 
          y.star.b <- y.star.b*(1-drop.mlt)
        }
        data.frame("y.new"=as.numeric(y.star.b), "Batch"=rep(names(n.batch)[bb], config.mat[bb, g]),
             "Group"=rep(names(n.group)[g], config.mat[bb, g]))
      })
      do.call("rbind", y.star.g)
    })
    Y.star <- do.call("rbind", Y.star)
    #if(any(is.na(Y.star))){print(i)}
  } 
  return(Y.star)
}


# Make y.star.b
#
# @param u copulas
# @param gg1 density object
#
# @return y.star.b
makeYstarB =  function(u, gg1){
  #Moved these steps outside the loop
  if(nrow(gg1)>1){
    difs<- diff(gg1$s)
    eps <- abs(mean(difs[is.finite(difs)], na.rm = TRUE))
    sapply(u, function(uu){
      yy  <- gg1$s[which.min(abs(gg1$Gy-uu))]
      yy  <- suppressWarnings(runif(1, yy-eps/2, yy+eps/2)) 
      yy
    })
  }else{
    numeric(length(u))
  } 
}