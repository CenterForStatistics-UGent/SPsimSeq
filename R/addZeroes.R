

addZeroes = function(Y, fracZero.logit.list, LL){
  if(mean(Y==min.val)>0.25){
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
}