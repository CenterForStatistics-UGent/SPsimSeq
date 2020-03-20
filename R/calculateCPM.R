# Calculates counts per millions of reads
calCPM <- function(X, const.mult = 1e6, norm.lib.size = TRUE, norm.factors = NULL, 
                   logt = FALSE, log_base = 2, prior.count = 1, ...){
  if(norm.lib.size){
    norm.factors = if(is.null(norm.factors)){
       edgeR::calcNormFactors(X)*colSums(X)
    } else{
     rep.int(1L, ncol(X))
    }
  }
  cpm <- X %*% diag(1/(norm.factors)) * const.mult
  if(logt){
    cpm <- log(cpm+prior.count, base=log_base)
  }
  return(cpm)
}