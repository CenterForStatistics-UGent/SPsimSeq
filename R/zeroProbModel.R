#' Predict zero probability using logistic rgression
#' @param cpm.data  log CPM matrix
#' @param logL log library size of the source data
#' @param zeroMat the matrix of zero indicators
#' @param n.mean.class a fraction of the number of genes
#' for the number of groups to be created for the mean log CPM of genes
#'
#' @return The coefficients of the estimated logistic regression
#' @importFrom stats glm.fit, binomial
#' @importFrom Hmisc cut2
zeroProbModel <- function(cpm.data, logL, zeroMat, n.mean.class){
  # Fit a logistic model for the probability of zeros, if there are enough of them
  mean.log.cpm <- rowMeans(cpm.data) #Means
  n.mean.class <- round(n.mean.class*nrow(cpm.data))
  mid.val <- as.numeric(as.character(cut2(mean.log.cpm, g = n.mean.class, 
                                          levels.mean = TRUE)))
  x1 <- rep(mid.val, times = length(logL))
  x2 <- rep(logL, each = length(mid.val))
  desMat = model.matrix(~x1*x2)
  fracZero.logit <- suppressWarnings(glm.fit(y = c(zeroMat), x = desMat, 
                                             family = binomial()))
  return(fracZero.logit$coef)
}