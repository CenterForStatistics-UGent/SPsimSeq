## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)   

library(SPsimSeq)

## ----load-packages, eval=FALSE, warning=FALSE, message=FALSE, echo=TRUE-------
#  # first install dependencies
#  pkg.depnd  <- c("SingleCellExperiment", "fitdistrplus",
#                 "edgeR", "Hmisc", "WGCNA", "limma", "mvtnorm", "plyr")
#  pkg.source <- c(1,2,1,2,2,1,2,2)
#  sapply(1:length(pkg.depnd), function(i){
#    pkg_available <- pkg.depnd[i] %in% rownames(installed.packages())
#    if(!pkg_available){
#      if(pkg.source[i] ==1){
#        BiocManager::install(pkg.depnd[i])
#      }
#      else{
#        install.packages(pkg.depnd[i])
#      }
#    }
#  })
#  
#  ## then install SPsimSeq
#  remotes::install_github("CenterForStatistics-UGent/SPsimSeq")
#  
#  
#  # load package
#  library(SPsimSeq)

## ---- eval=TRUE, warning=FALSE, fig.width=8, fig.height=4---------------------
 # load the Zhang bulk RNA-seq data (availabl with the package) 
 data("zhang.data.sub") 
 
 # filter genes with sufficient expression (important step to avoid bugs) 
 zhang.counts <- zhang.data.sub$counts[rowSums(zhang.data.sub$counts > 0)>=5, ]  
 MYCN.status  <- zhang.data.sub$MYCN.status 
 
 set.seed(6452)
 zhang.counts2 <- zhang.counts[sample(nrow(zhang.counts), 3000), ]
 # We simulate only a single data (n.sim = 1) with the following property
 # - 3000 genes ( n.genes = 3000) 
 # - 172 samples (tot.samples = 172) -- equal to the source data
 # - the samples are equally divided into 2 groups each with 90 samples 
 #   (group.config = c(0.5, 0.5)) -- almost equal to the source data
 # - all samples are from a single batch (batch = NULL, batch.config = 1)
 # - we add 10% DE genes (pDE = 0.1) 
 # - the DE genes have a log-fold-change of at least 0.5 in 
 #   the source data (lfc.thrld = 0.5)
 # - we do not model the zeroes separately, they are the part of density 
 #    estimation (model.zero.prob = FALSE)
 
 # simulate data
 sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = zhang.counts2, batch = NULL,
                          group = MYCN.status, n.genes = 3000, batch.config = 1,
                          group.config = c(0.5, 0.5), tot.samples = 172, 
                          pDE = 0.1, lfc.thrld = 0.5, result.format = "list",  
                          seed = 2581988)
 sim.data.bulk1 <- sim.data.bulk[[1]]
 head(sim.data.bulk1$counts[, 1:5])  # count data
 head(sim.data.bulk1$colData)        # sample info
 head(sim.data.bulk1$rowData)        # gene info
 
 
# compare the distributions of the mean expressions, variability, 
# and fraction of zero counts per gene

 library(LSD) # for generating heatmap plots
 
# normalize counts for comparison  
Y0.log.cpm <- log2(edgeR::cpm(zhang.counts2)+1)
Y1.log.cpm <- log2(edgeR::cpm(sim.data.bulk1$counts)+1)
Y0.log.cpm <- Y0.log.cpm[rowMeans(Y0.log.cpm>0)>=0.1, ]
Y1.log.cpm <- Y1.log.cpm[rowMeans(Y1.log.cpm>0)>=0.1, ]
 

rowVars <- function(X){apply(X, 1, var, na.rm=TRUE)}
rowCVs <- function(X){apply(X, 1, function(x) sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))}

 
par(mfrow=c(1, 3))
boxplot(list(real.data=log(colSums(zhang.counts2)), 
             simulated.data=log(sim.data.bulk1$colData$sim.Lib.Size)), 
        main="library size") 
boxplot(list(real.data=rowMeans(Y0.log.cpm), 
             simulated.data=rowMeans(Y1.log.cpm)), 
        main="mean expression of genes") 
boxplot(list(real.data=rowVars(Y0.log.cpm), 
             simulated.data=rowVars(Y1.log.cpm)), 
        main="variance of gene expressions") 
 
 

# compare the relationship between the mean and variability
par(mfrow=c(1,3), mar=c(4,4,4,1))
heatscatter(rowMeans(Y0.log.cpm), rowCVs(Y0.log.cpm), ylim=c(0, 6), xlim=c(0, 16),
            colpal="bl2gr2rd", main="real data", xlab="mean log2-CPM", 
            ylab="coefficients of variation", cexplot=0.5, alpha = 60, cex.lab=1.25)
heatscatter(rowMeans(Y1.log.cpm), rowCVs(Y1.log.cpm), ylim=c(0, 6), xlim=c(0, 16),
     main="SPsimSeq", xlab="mean log2-CPM", ylab="coefficients of variation", 
     cexplot=0.5, alpha = 60, colpal="bl2gr2rd", cex.lab=1.25)


n.gride <- 1000
min.g   <- seq(0, 20, length.out = n.gride+1)[-n.gride]
max.g   <- seq(0, 20, length.out = n.gride+1)[-1] 
mid.g   <- (min.g+max.g)/2
f.real  <- sapply(1:n.gride, function(r){
  x <- Y0.log.cpm[rowMeans(Y0.log.cpm)<=max.g[r] & rowMeans(Y0.log.cpm)>min.g[r],]
  y <- ifelse(!is.null(dim(x)), mean(rowCVs(x)), mean(sd(x)/mean(x))) 
  y
})
f.SPsim <- sapply(1:n.gride, function(r){
  x <- Y1.log.cpm[rowMeans(Y1.log.cpm)<=max.g[r] & rowMeans(Y1.log.cpm)>min.g[r],]
  y <- ifelse(!is.null(dim(x)), mean(rowCVs(x)), mean(sd(x)/mean(x))) 
  y
})


sm1 <- loess(I(f.SPsim-f.real)~mid.g) 
plot(mid.g, f.SPsim-f.real, xlim=c(0, 14), col="lightskyblue", pch=20, cex.lab=1.25,
     cex.main=1.4, main="SPsimSeq - real data", ylab="difference", xlab="mean log2-CPM")
lines(mid.g,predict(sm1, newdata = mid.g), col="blue", lwd=3) 


# compare the correlation between genes and samples 
cor.mat.Y0 <- cor(t(Y0.log.cpm))
cor.mat.Y1 <- cor(t(Y1.log.cpm)) 

cor.vec.Y0 <- cor.mat.Y0[upper.tri(cor.mat.Y0)]
cor.vec.Y1 <- cor.mat.Y1[upper.tri(cor.mat.Y1)] 

par(mfrow=c(1,3), mar=c(4,4,3.5,1))
hist(cor.vec.Y0, nclass = 30, probability = TRUE, 
     border="gray", col="steelblue1", main="real data", xlab="pairwise correlation between genes", 
     ylim=c(0, 3.5), xlim=c(-1, 1), cex.lab=1.25)
hist(cor.vec.Y1, nclass = 30, probability = TRUE, border="gray",
     col="steelblue1",  main="SPsimSeq", xlab="pairwise correlation between genes",
     ylim=c(0, 3.5), xlim=c(-1, 1), cex.lab=1.25)

plot(seq(-1, 1, 0.1), seq(-1, 1, 0.1), type="n", xlab="quantile (real data)", 
     ylab="quantile (simulated data)",  main="correlation quantile-quantile plot")
abline(0, 1, col="gray")
points(quantile(cor.vec.Y0, seq(0, 1, 0.001)), quantile(cor.vec.Y1, seq(0, 1, 0.001)), 
       col="blue", pch=20, cex=1.5, cex.lab=1.25)  

