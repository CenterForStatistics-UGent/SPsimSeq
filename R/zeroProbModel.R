#' predict zero probability using logistic rgression
#' @param cpm.data  log CPM matrix
#' @param L library size of the source data
#' @param n.mean.class a fraction of the number of genes
#' for the number of groups to be created for the mean log CPM of genes
#' @param subset.data a logical value to subset columns(samples) if the size of the data is too big
#' for space and computation time saving
#' @param n.samples an integer indicating the number of samples/cells to sample (if subset.data=TRUE)
#' @param const a small constant (>0) to be added to the CPM before log transformation, to avoid  log(0).
# default is 1e-5
#' @param  ... further arguments passed to or from other methods. 
#'
#' @return a GLM class object contating the estimated logistic regression
#' @importFrom stats glm
#' @importFrom Hmisc cut2
zeroProbModel <- function(cpm.data, L, n.mean.class=0.2, subset.data=FALSE, 
                          n.samples=400, const){
  
  if(subset.data & (ncol(cpm.data)>n.samples)){
    sel.cols <- sample(ncol(cpm.data), n.samples)
    cpm.data <- cpm.data[, sel.cols]
    L <- L[sel.cols]
  } 
  
  # fit logistic model for the probability of zeros
  mean.log.cpm <- rowMeans(cpm.data)
  y <- c(as.integer(cpm.data==log(const)))
  n.mean.class <- round(n.mean.class*nrow(cpm.data))
  mid.val <- as.numeric(as.character(cut2(mean.log.cpm, g=n.mean.class, levels.mean = TRUE)))
  log.L <- log(L)
  x1 <- rep(mid.val, times=length(log.L))
  x2 <- rep(log.L, each=length(mid.val))
  fracZero.logit <- suppressWarnings(glm(y~x1*x2, family = "binomial"))
  fracZero.logit
}