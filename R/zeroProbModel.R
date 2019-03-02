#' predict zero probability using logistic rgression
#'
#' More detailed description
#'
#' @param count.data a count data matrix with rows represent genes and columns represent cells/samples
#' @param n.mean.class an integer for the number of groups to be created for the mean log CPM of genes
#' @param subset.data a logical value to subset columns(samples) if the size of the data is too big
#' for space and computation time saving
#' @param n.samples an integer indicating the number of samples/cells to sample (if subset.data=TRUE)
#' @param  ... further arguments passed to or from other methods. 
#'
#' @return a GLM class object contating the estimated logistic regression
#' @examples
#'  # example
#'  #count.data <- 
#'  #fracZero.logit <- zeroProbModel()
#'  #summary(fracZero.logit)
#'  #par(mfrow=c(1,2))
#'  #plot(exp(log.L), predict(fracZero.logit, type="response",                        
#'  #newdata=data.frame(x1=median(mean.log.cpm), x2=log.L)), ylim=c(0,1))
#'  #plot(mean.log.cpm, predict(fracZero.logit, type="response",
#'  #newdata=data.frame(x1=mean.log.cpm, x2=median(log.L))), ylim=c(0,1))
#'  
#' @export 
zeroProbModel <- function(count.data, n.mean.class=30, subset.data=FALSE, n.samples=400, ...){

  if(subset.data & (ncol(count.data)>n.samples)){
    count.data <- count.data[, sample(ncol(count.data), n.samples)]
  }
  count.data <- count.data[rowSums(count.data>0)>5, ]
  cpm.data <- log(calCPM(count.data)+1)

  # fit logistic model for the probability of zeros
  mean.log.cpm <-rowMeans(cpm.data)
  Z <- ifelse(cpm.data==0, 1, 0)
  mid.val <- as.numeric(as.character(cut2(mean.log.cpm, g=n.mean.class, levels.mean = TRUE)))
  # mid.val <- hist(mean.log.cpm, nclass = n.mean.class, plot = FALSE)
  # breaks <- mid.val$breaks
  # mid.val <- as.numeric(as.character(cut2(mean.log.cpm, cuts = breaks, levels.mean = TRUE)))
  log.L <- log(colSums(count.data))
  y  <- as.vector(Z)
  x1 <- rep(mid.val, times=length(log.L))
  x2 <- rep(log.L, each=length(mid.val))
  #fracZero.logit <- glm(y~I(x1)+I(x2)+I(x2^2)+I(x1*x2)+I(x1*x2^2), family = "binomial")
  fracZero.logit <- glm(y~x1*x2, family = "binomial")
  fracZero.logit
}

# count.data <- counts(s.data) 
# summary(fracZero.logit)
# par(mfrow=c(1,2))
# plot(exp(log.L), predict(fracZero.logit, type="response",                        
#                          newdata=data.frame(x1=median(mean.log.cpm), x2=log.L)), ylim=c(0,1))
# plot(mean.log.cpm, predict(fracZero.logit, type="response",
#                            newdata=data.frame(x1=mean.log.cpm, x2=median(log.L))), ylim=c(0,1))
# 

