#' Predict zero probability using logistic rgression
#' @param cpm.data  log CPM matrix
#' @param logL log library size of the source data
#' @param n.mean.class a fraction of the number of genes
#' for the number of groups to be created for the mean log CPM of genes
#' @param subset.data a logical value to subset columns(samples) if the size of the data is too big
#' for space and computation time saving
#' @param n.samples an integer indicating the number of samples/cells to sample (if subset.data=TRUE)
#' @param const a small constant (>0) to be added to the CPM before log transformation, to avoid  log(0)
#'
#' @return The coefficients of the estimated logistic regression
#' @importFrom stats glm.fit, binomial
#' @importFrom Hmisc cut2
zeroProbModel <- function(cpm.data, logL, zeroMat, n.mean.class = 0.2, subset.data = FALSE, 
                          n.samples = 400){
  #Subset data if needed
  if(subset.data & (ncol(cpm.data) > n.samples)){
    sel.cols <- sample(ncol(cpm.data), n.samples)
    cpm.data <- cpm.data[, sel.cols]
    logL <- logL[sel.cols]
    zeroMat = zeroMat[, sel.cols]
  } 
  # Fit a logistic model for the probability of zeros
  mean.log.cpm <- rowMeans(cpm.data) #Means
  n.mean.class <- round(n.mean.class*nrow(cpm.data))
  mid.val <- as.numeric(as.character(cut2(mean.log.cpm, g = n.mean.class, 
                                          levels.mean = TRUE)))
  x1 <- rep(mid.val, times = length(log.L))
  x2 <- rep(log.L, each = length(mid.val))
  desMat = model.matrix(~x1*x2)
  fracZero.logit <- suppressWarnings(glm.fit(y = c(zeroMat), x = desMat, family = binomial()))
  fracZero.logit$coef
}