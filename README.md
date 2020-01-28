---
title: "SPsimSeq (v2.0.1)"
output: "github_document"
---

This is the github repo for the SPsimSeq R package. 

# Overview
SPsimSeq uses a specially designed exponential family for density estimation to constructs the distribution of gene expression levels from a given real RNA sequencing data (single-cell or bulk), and subsequently, simulates a new dataset from the estimated marginal distributions using Gaussian-copulas to retain the dependence between genes. It allows simulation of multiple groups and batches with any required sample size and library size.




 
# Installation

```
### first install dependencies 
pkg.depnd  <- c("SingleCellExperiment", "fitdistrplus",  
               "edgeR", "Hmisc", "WGCNA", "limma", "mvtnorm")
pkg.source <- c(1,2,1,2,2,1,2)
sapply(1:length(pkg.depnd), function(i){
  pkg_available <- pkg.depnd[i] %in% rownames(installed.packages())
  if(!pkg_available){ 
    if(pkg.source[i] ==1){
      BiocManager::install(pkg.depnd[i])
    }
    else{
      install.packages(pkg.depnd[i])
    } 
  }
}) 


### then install SPsimSeq
remotes::install_github("CenterForStatistics-UGent/SPsimSeq")


### load package
library(SPsimSeq) 
```
 
# Demonstrations
 
## Example 1: simulating bulk RNA-seq

``` 
# load the Zhang bulk RNA-seq data (availabl with the package) 
data("zhang.data.sub") 

# filter genes with sufficient expression (important step) 
zhang.counts <- zhang.data.sub$counts[rowSums(zhang.data.sub$counts > 0)>=5, ]  
MYCN.status  <- zhang.data.sub$MYCN.status 


# We simulate only a single data (n.sim = 1) with the following property
# - 2000 genes ( n.genes = 2000) 
# - 20 samples (tot.samples = 20) 
# - the samples are equally divided into 2 groups each with 90 samples 
#   (group.config = c(0.5, 0.5))
# - all samples are from a single batch (batch = NULL, batch.config = 1)
# - we add 10% DE genes (pDE = 0.1) 
# - the DE genes have a log-fold-change of at least 0.5 in 
#   the source data (lfc.thrld = 0.5)
# - we do not model the zeroes separately, they are the part of density 
#    estimation (model.zero.prob = FALSE)

# simulate data
set.seed(6452)
zhang.counts2 <- zhang.counts[sample(nrow(zhang.counts), 2000), ]
sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = zhang.counts2, batch = NULL,
                          group = MYCN.status, n.genes = 2000, batch.config = 1,
                          group.config = c(0.5, 0.5), tot.samples = 20, 
                          pDE = 0.1, lfc.thrld = 0.5, 
                          result.format = "list",  seed = 2581988)

sim.data.bulk1 <- sim.data.bulk[[1]]                              
head(sim.data.bulk1$counts[, seq_len(5)])  # count data
head(sim.data.bulk1$colData)               # sample info
head(sim.data.bulk1$rowData)               # gene info
```


## Example 2: simulating single cell RNA-seq from a single batch (read-counts)

```
# we simulate only a single scRNA-seq data (n.sim = 1) with the following property
# - 2000 genes (n.genes = 2000) 
# - 100 cells (tot.samples = 100) 
# - the cells are equally divided into 2 groups each with 50 cells 
#   (group.config = c(0.5, 0.5))
# - all cells are from a single batch (batch = NULL, batch.config = 1)
# - we add 10% DE genes (pDE = 0.1) 
# - the DE genes have a log-fold-change of at least 0.5
# - we model the zeroes separately (model.zero.prob = TRUE)
# - the ouput will be in SingleCellExperiment class object (result.format = "SCE")


library(SingleCellExperiment)

# load the NGP nutlin data (availabl with the package, processed with 
# SMARTer/C1 protocol, and contains read-counts)
data("scNGP.data")

# filter genes with sufficient expression level (important step) 
scNGP.data2 <- scNGP.data[rowSums(counts(scNGP.data) > 0)>=5, ]  
treatment <- ifelse(scNGP.data2$characteristics..treatment=="nutlin",2,1)

set.seed(654321)
scNGP.data2 <- scNGP.data2[sample(nrow(scNGP.data2), 2000), ]

# simulate data (we simulate here only a single data, n.sim = 1)
sim.data.sc <- SPsimSeq(n.sim = 1, s.data = scNGP.data2, batch = NULL,
                        group = treatment, n.genes = 2000, batch.config = 1,
                        group.config = c(0.5, 0.5), tot.samples = 100, 
                        pDE = 0.1, lfc.thrld = 0.5, model.zero.prob = TRUE,
                        result.format = "SCE",  seed = 2581988)

sim.data.sc1 <- sim.data.sc[[1]]
class(sim.data.sc1)
head(counts(sim.data.sc1)[, seq_len(5)])
colData(sim.data.sc1)
rowData(sim.data.sc1)
```

