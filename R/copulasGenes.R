
#Estimates gausian copulas
#
#@param X cpm.data
#@param batch batch
#@param n.batch n.batch
#@return
#
#
#@importFrom mvtnorm rmvnorm
#@importFrom WGCNA cor
#@importFrom limma voom
genesCopula <- function(X, batch, n.batch){
  cpl.list <- lapply(sort(unique(batch)), function(bb){ 
    voom.res <- limma::voom(X[, batch==bb])
    W <- voom.res$weights
    cor.mat <- WGCNA::cor(x=t(X[, batch==bb]), weights.x = t(W), 
                          use = "pairwise.complete.obs")
    Z       <- mvtnorm::rmvnorm(n=n.batch[bb], sigma = cor.mat)  
    Cpl     <- apply(Z, 2, function(z) pnorm(z))
    colnames(Cpl) <- rownames(X)
    Cpl
  }) 
  cpl.list
}

