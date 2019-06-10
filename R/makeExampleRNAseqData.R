#' Make an example RNA-seq data.
#' 
#' @param n.gene number of genes.
#' @param n.sample number of samples/cells.
#' @param n.group number of biological group.
#' @param n.batch number of batchs.
#' @param p.DEgenes fraction of genes with DE between simulated groups
#' @param common.BCV common biological coefficient of variation
#' @param type type gene expression data, either bulk RNA-seq (type="bulk") or 
#' single cell RNA-seq (type="sc")
#' @return a SingleCellExperiment object.
#'  
#' @examples
#' # make.example.data
#' @export
#' @importFrom stats rlnorm rnbinom 
#' @importFrom SingleCellExperiment counts colData rowData SingleCellExperiment

make.example.data <- function(n.gene=2000, n.sample=30, n.group=2, n.batch=3, p.DEgenes=0.2,
                              common.BCV=0.25, type="bulk"){ 
  # relative abundance
  r0 <- rgamma(n.gene, 5, 1)/sum(rgamma(n.gene, 5, 1))
  
  # add biological signal
  group.fact <- matrix(1, ncol=n.group, nrow=n.gene)
  for(i in 1:ncol(group.fact)){
    x <- rlnorm(p.DEgenes*n.gene/2, meanlog = 0.5, sdlog = 0.4)
    group.fact[sample(nrow(group.fact), round(p.DEgenes*n.gene/2)),i] <- x
  } 
  
  # add batch effect in each biological group
  batch.fact <- matrix(1, ncol=n.batch, nrow=n.gene)
  for(i in 1:ncol(batch.fact)){
    x <- rlnorm(n.gene/2, meanlog = 0.5, sdlog = 0.4)
    batch.fact[sample(nrow(batch.fact), round(n.gene/2)),i] <- x
  }
  
  r1 <- do.call(cbind, lapply(1:ncol(group.fact), function(i){
    do.call(cbind, lapply(1:ncol(batch.fact), function(j){ 
      replicate(n.sample, r0*group.fact[, i]*batch.fact[, j])
    })) 
  }))
  
  LS <- rlnorm(n.sample*n.group*n.batch, meanlog = log(25*n.gene), sdlog = 0.1)
  group <- gl(n.group, n.sample*n.batch)
  batch <- rep(gl(n.batch, n.sample), times=n.group)
  
  mu.dat <- r1%*%diag(LS)
  disp   <- (common.BCV+ 1/sqrt(rowMeans(mu.dat))) #rchisq(n.gene, 2)
  cnt <- t(sapply(1:nrow(mu.dat), function(gn){
    sapply(1:ncol(mu.dat), function(cl){
      rnbinom(1, mu=mu.dat[gn, cl], size=1/disp[gn])
    })
  }))
  rownames(cnt) <- paste0("Gene_", 1:n.gene)
  colnames(cnt) <- paste0("Sample_", 1:ncol(cnt))
  col.dat <- data.frame(E.LS=LS, Group=group, Batch=batch, 
                        row.names = colnames(cnt))
  
  group.fact <- as.data.frame(group.fact)
  colnames(group.fact) <- paste0("DE.fac", 1:n.group)
  
  batch.fact <- as.data.frame(batch.fact)
  colnames(batch.fact) <- paste0("Batch.fac", 1:n.batch)
  
  row.dat <- data.frame(group.fact, batch.fact, rho=r0, row.names = rownames(cnt))
  
  
  if(type=="bulk"){
    SingleCellExperiment(assays = list(counts=cnt), 
                         colData = col.dat, rowData=row.dat) 
  }
  else if(type=="sc"){
    # # add excess zero
    # ## model zero prob
    # logCPM <- log(calCPM(cnt)+1)
    # Z = as.vector(logCPM<quantile(logCPM, 0.25))
    # logLS = log(colSums(cnt))
    # meanLogCPM <- rowMeans(logCPM)
    # 
    # x1 <- rep(meanLogCPM, times=length(logLS))
    # x2 <- rep(logLS, each=length(meanLogCPM))
    # 
    # lr.fit <- glm(Z~x1*x2, family = "binomial")
    SingleCellExperiment(assays = list(counts=cnt), 
                         colData = col.dat, rowData=row.dat) 
  }
}


