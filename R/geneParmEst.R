# Gene level param estimates for density estimation
gene.parm.est <- function(cpm.data.i, batch, group, null.group,
                          sub.batchs, de.ind, model.zero.prob,
                          min.val, w,...){
  batch.est <- lapply(sub.batchs, function(b){ 
    if(de.ind==0){
      Y0 <- cpm.data.i[(batch==b & group==null.group)] 
      Y = if(model.zero.prob & mean(Y0==min.val)>0.25){
        Y0[Y0>min.val]
        } else Y0
      
    if(sum(Y>min.val)>3){
        countY   <- obtCount(Y=Y, w=w)
        parm.est <- fitLLmodel(yy=countY)
        # if(!is.null(parm.est$betas) & length(countY$S)>=3){
        #   parm.est
        # } 
        # else{NULL}  
        parm.est
      } 
      else{
        list(S=0, lls=0, uls=0, Ny=0, Y=Y, mu.hat=mean(Y), sig.hat=sd(Y), 
                w=w, betas=parmEstOut(llm = NULL)$beta.hat.vec, 
                v = parmEstOut(llm = NULL)$V.hat.mat)
        }
    }
    else{
      Y0 <- split(cpm.data.i[batch==b], group[batch==b])
      Y <- lapply(Y0, function(y0){
        if(model.zero.prob & mean(y0==min.val)>0.25) {y0[y0>min.val]}
        else {y0}
      }) 
      lapply(Y, function(y){
        if(sum(y>min.val)>3){
          countY   <- obtCount(Y=y, w=w)
          parm.est <- fitLLmodel(yy=countY)
          # if(!is.null(parm.est$betas) & length(countY$S)>=3){
          #   parm.est
          # } 
          # else{NULL}  
          parm.est
        } 
        else{
          list(S=0, lls=0, uls=0, Ny=0, Y=y, mu.hat=mean(y), sig.hat=sd(y), 
               w=w, betas=parmEstOut(llm = NULL)$beta.hat.vec, 
               v = parmEstOut(llm = NULL)$V.hat.mat)
        }
      }) 
    }
  }) 
  return(batch.est)
}