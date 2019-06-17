# predict zero probability using logistic rgression
#
# More detailed description
#
# @param cpm.data  log CPM matrix
# @param L library size of the source data
# @param n.mean.class a fraction of the number of genes
# for the number of groups to be created for the mean log CPM of genes
# @param subset.data a logical value to subset columns(samples) if the size of the data is too big
# for space and computation time saving
# @param n.samples an integer indicating the number of samples/cells to sample (if subset.data=TRUE)
# @param simCtr (integer) seed number. Used if subset.data=TRUE
# @param const a small constant (>0) to be added to the CPM before log transformation, to avoid  log(0).
# default is 1e-5
# @param  ... further arguments passed to or from other methods. 
#
# @return a GLM class object contating the estimated logistic regression
# @examples
#  # example
#  #count.data <- 
#  #fracZero.logit <- zeroProbModel()
#  #summary(fracZero.logit)
#  #par(mfrow=c(1,2))
#  #plot(exp(log.L), predict(fracZero.logit, type="response",                        
#  #newdata=data.frame(x1=median(mean.log.cpm), x2=log.L)), ylim=c(0,1))
#  #plot(mean.log.cpm, predict(fracZero.logit, type="response",
#  #newdata=data.frame(x1=mean.log.cpm, x2=median(log.L))), ylim=c(0,1))
#  
# @export  
# @importFrom stats glm
# @importFrom Hmisc cut2
zeroProbModel <- function(cpm.data, L, n.mean.class=0.2, subset.data=FALSE, 
                          n.samples=400, simCtr, const, ...){
  
  if(subset.data & (ncol(cpm.data)>n.samples)){
    #set.seed(simCtr)
    sel.cols <- sample(ncol(cpm.data), n.samples)
    cpm.data <- cpm.data[, sel.cols]
    L <- L[sel.cols]
    #set.seed(NULL)
  } 
  
  # fit logistic model for the probability of zeros
  mean.log.cpm <- rowMeans(cpm.data)
  Z <- ifelse(cpm.data==log(const), 1, 0)
  n.mean.class <- round(n.mean.class*nrow(cpm.data))
  mid.val <- as.numeric(as.character(cut2(mean.log.cpm, g=n.mean.class, levels.mean = TRUE)))
  # mid.val <- hist(mean.log.cpm, nclass = n.mean.class, plot = FALSE)
  # breaks <- mid.val$breaks
  # mid.val <- as.numeric(as.character(cut2(mean.log.cpm, cuts = breaks, levels.mean = TRUE)))
  log.L <- log(L)
  y  <- as.vector(Z)
  x1 <- rep(mid.val, times=length(log.L))
  x2 <- rep(log.L, each=length(mid.val))
  #fracZero.logit <- glm(y~I(x1)+I(x2)+I(x2^2)+I(x1*x2)+I(x1*x2^2), family = "binomial")
  fracZero.logit <- suppressWarnings(glm(y~x1*x2, family = "binomial"))
  fracZero.logit
}
 

# summary(fracZero.logit)
# par(mfrow=c(1,2))
# plot(exp(log.L), predict(fracZero.logit, type="response",
#                          newdata=data.frame(x1=max(mean.log.cpm), x2=log.L)))
# plot(mean.log.cpm, predict(fracZero.logit, type="response",
#                            newdata=data.frame(x1=mean.log.cpm, x2=min(log.L))), ylim=c(0,1))


