# A function to obtain copulas or uniform random variables
obtCopulasBatch <- function(genewiseCor, cpm.data, batch, n.batch){
  if(genewiseCor){
    copulas.batch <- genesCopula(X = cpm.data, batch = batch, n.batch=n.batch)
  }
  else{
    copulas.batch <- lapply(sort(unique(batch)), function(bb){ 
      U    <- matrix(runif(n.batch[bb]*nrow(cpm.data), 0, 1), n.batch[bb], nrow(cpm.data))   
      colnames(U) <- rownames(cpm.data)
      U
    }) 
  }
  copulas.batch
}