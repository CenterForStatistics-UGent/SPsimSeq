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
    zeroModel = zeroProbModel(cpm.data = cpm.data[, coln], logL = log(LS[coln]), 
                  zeroMat[, coln],  n.mean.class = n.mean.class, 
                  minFracZeroes = minFracZeroes)
    #Calculate zero fractions of each gene within the batches
    zeroFrac = rowMeans(zeroMat[, coln])
    #Calculate gene-wise means
    geneMeans = rowMeans(cpm.data[, coln])
    #Retain the means of the genes exceeding the zero fraction threshold
    meansLarge = geneMeans[zeroFrac > minFracZeroes]
    list(zeroModel = zeroModel, meansLarge = meansLarge)
  })
}