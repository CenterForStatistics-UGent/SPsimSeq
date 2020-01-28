# A wrapper function to call the ZerpProbModel
fracZeroLogitModel <- function(s.data, cpm.data, batch, const, sub.batchs, ...){
  if(is(s.data, "SingleCellExperiment")){
    fracZero.logit.list <- lapply(unique(sub.batchs), function(b){
      zeroProbModel(cpm.data = cpm.data[, batch==b], L=colSums(counts(s.data)[, batch==b]),
                    simCtr=NULL, const=const, ...)
    }) 
  }
  else if(is(s.data, "data.frame") | is(s.data, "matrix")){
    fracZero.logit.list <- lapply(unique(sub.batchs), function(b){
      zeroProbModel(cpm.data = cpm.data[, batch==b], L=colSums(s.data[, batch==b]), 
                    simCtr=NULL, const=const,...)
    })
  }
  fracZero.logit.list
}