# A wrapper function to call the ZerpProbModel
fracZeroLogitModel <- function(s.data, cpm.data, batch, const, sub.batchs){
    lapply(unique(sub.batchs), function(b){
      zeroProbModel(cpm.data = cpm.data[, batch==b], L=colSums(s.data[, batch==b]), 
                   const=const)
    })
}