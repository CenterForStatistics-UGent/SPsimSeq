#' Extract data and iterate over batches to estimate zero probability models
#'
#' @param s.data,cpm.data raw and transformed data 
#' @param batch the batch vector
#' @param n.mean.class see zeroProbModel
#' @param minFracZeroes minimum fraction of zeroes before zero-inflation is applied
#'
#' @return a list of binomial regression parameters
fracZeroLogitModel <- function(s.data, batch, cpm.data, n.mean.class, 
                               minFracZeroes){
  LS = colSums(s.data)
  zeroMat = s.data == 0
  zeroModels = tapply(colnames(s.data), batch, function(coln){
    zeroProbModel(cpm.data = cpm.data[, coln], logL = log(LS[coln]), 
                  zeroMat[, coln],  n.mean.class = n.mean.class, 
                  minFracZeroes = minFracZeroes)
  })
  zeroFracs = vapply(colnames(s.data), function(gene){
    tapply(zeromat[gene, ], batch, mean)
  }, FUN.VALUE = numeric(length(unique(batch))))
  list(zeroModels = zeroModels, exceedIds = exceedIds)
}