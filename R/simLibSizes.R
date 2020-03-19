#' A function simulate library sizes 
#'
#' @param sub.batchs 
#' @param LS observed library sizes
#' @param lib.size.params Parameters of the lognormal
#' @param n.batch number of batches
#' @param batch batches
#' @param n.group number of groups
#' @param config.mat groups
#' 
#' @return simulated library sizes
simLibSize <- function(sub.batchs, LS, lib.size.params, n.batch, batch, n.group, config.mat){
  LL <- lapply(seq_along(sub.batchs), function(b){
    if(is.null(lib.size.params)){
      L.b      <- LS[batch==sub.batchs[b]]
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