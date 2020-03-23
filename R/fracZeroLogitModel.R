#' Extract data and iterate over batches to estimate zero probability models
#'
#' @param s.data 
#' @param batch 
#'
#' @return a list of binomial regression parameters
fracZeroLogitModel <- function(s.data, batch, cpm.data){
  LS = colSums(s.data)
  zeroMat = s.data == 0
  tapply(colnames(s.data), batch, function(coln){
    zeroProbModel(cpm.data = cpm.data[, coln], logL = log(LS[coln]), 
                  zeroMat[, coln])
  })
}