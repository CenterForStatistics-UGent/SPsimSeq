# A function simulate library size 
simLibSize <- function(sub.batchs, L, lib.size.params, n.batch, batch, n.group, config.mat, ...){
  LL <- lapply(seq_len(length(sub.batchs)), function(b){
    if(is.null(lib.size.params)){
      L.b      <- L[batch==sub.batchs[b]]
      fit.ln   <- fitdist(as.numeric(L.b), distr = "lnorm")$estimate 
      L.b.pred <- rlnorm(n.batch[b], fit.ln[["meanlog"]], fit.ln[["sdlog"]])
    }
    else{
      if(length(lib.size.params) != 2 | is.null(names(lib.size.params))){
        stop("The log-normal parameters for the distribution of library sizes must be submitted in a named vector of size 2. 
             Example, lib.size.params = c(meanlog=10, sdlog=0.2). See also ?rlnorm()")
      }else{
        L.b.pred <- rlnorm(n.batch[b], lib.size.params[["meanlog"]], lib.size.params[["sdlog"]])
      } 
    }
    
    
    ## randomly split into the groups
    gr <- rep(seq_len(length(n.group)), config.mat[b, ])
    split(L.b.pred, gr) 
  })
  return(LL)
}