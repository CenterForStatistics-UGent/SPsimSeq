SPsimSeq
================

This is the github repo for the SPsimSeq R package.

Overview
========

SPsimSeq uses a specially designed exponential family for density estimation to constructs the distribution of gene expression levels from a given real RNA sequencing data (single-cell or bulk), and subsequently, simulates a new dataset from the estimated marginal distributions using Gaussian-copulas to retain the dependence between genes. It allows simulation of multiple groups and batches with any required sample size and library size.

Installation
============

Github installation

``` r
library(devtools)
install_github("CenterForStatistics-UGent/SPsimSeq")
```

BioConductor installation

``` r
library(BiocManager)
BiocManager::install("SPsimSeq", update = FALSE)
```

``` r
suppressPackageStartupMessages(library(SPsimSeq))
cat("SPsimSeq package version", 
    as.character(packageVersion("SPsimSeq")), "\n")
```

    ## SPsimSeq package version 0.99.0

Demonstrations
==============

Example 1: simulating bulk RNA-seq
----------------------------------

``` r
# load the Zhang bulk RNA-seq data
data("zhang.data.sub") 
# filter genes with sufficient expression (important step) 
zhang.counts <- zhang.data.sub$counts[rowSums(zhang.data.sub$counts > 0)>=5, ]  
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
sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = zhang.counts2,
                          group = zhang.data.sub$MYCN.status, n.genes = 2000, 
                          batch.config = 1,
                          group.config = c(0.5, 0.5), tot.samples = 20, 
                          pDE = 0.1, lfc.thrld = 0.5, 
                          result.format = "list")
```

    ## Estimating featurewise correlations ...

    ## 

    ## Selecting candidate DE genes ...

    ## Note: The number of DE genes detected in the source data is 79 and the number of DE genes required to be included in the simulated data is 200. Therefore, candidiate DE genes are sampled with replacement.

    ## Fitting zero probability model ...

    ## Estimating densities ...

    ## Configuring design ...

    ## Simulating data ...

    ##  ...1 of 1

``` r
sim.data.bulk1 <- sim.data.bulk[[1]]                              
head(sim.data.bulk1$counts[, seq_len(5)])  # count data
```

    ##        Sample_1 Sample_2 Sample_3 Sample_4 Sample_5
    ## Gene_1        0        0        0        0        7
    ## Gene_2        1        1        1        8        0
    ## Gene_3        2       38        0       10       28
    ## Gene_4        3       25        2      132      232
    ## Gene_5        0        0        8        7       12
    ## Gene_6        1        1        0       58        1

``` r
head(sim.data.bulk1$colData)               # sample info
```

    ##          Batch Group sim.Lib.Size
    ## Sample_1     1     1      1484907
    ## Sample_2     1     1      1252074
    ## Sample_3     1     1      1861492
    ## Sample_4     1     1      1644716
    ## Sample_5     1     1      1769902
    ## Sample_6     1     1      1260043

``` r
head(sim.data.bulk1$rowData)               # gene info
```

    ##        DE.ind     source.ID
    ## Gene_1   TRUE     HIST1H2BM
    ## Gene_2   TRUE RP11-378A13.1
    ## Gene_3   TRUE         HTR2C
    ## Gene_4   TRUE          ZIC5
    ## Gene_5   TRUE           WT1
    ## Gene_6   TRUE       C2orf81

Example 2: simulating single cell RNA-seq from a single batch (read-counts)
---------------------------------------------------------------------------

``` r
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
suppressPackageStartupMessages(library(SingleCellExperiment))
# load the NGP nutlin data (availabl with the package, processed with 
# SMARTer/C1 protocol, and contains read-counts)
data("scNGP.data")
# filter genes with sufficient expression level (important step) 
scNGP.data2 <- scNGP.data[rowSums(counts(scNGP.data) > 0)>=5, ]  
treatment <- ifelse(scNGP.data2$characteristics..treatment=="nutlin",2,1)
set.seed(654321)
scNGP.data2 <- scNGP.data2[sample(nrow(scNGP.data2), 2000), ]
# simulate data (we simulate here only a single data, n.sim = 1)
sim.data.sc <- SPsimSeq(n.sim = 1, s.data = scNGP.data2,
                        group = treatment, n.genes = 2000, batch.config = 1,
                        group.config = c(0.5, 0.5), tot.samples = 100, 
                        pDE = 0.1, lfc.thrld = 0.5, model.zero.prob = TRUE,
                        result.format = "SCE")
```

    ## Estimating featurewise correlations ...

    ## Selecting candidate DE genes ...

    ## Note: The number of DE genes detected in the source data is 51 and the number of DE genes required to be included in the simulated data is 200. Therefore, candidiate DE genes are sampled with replacement.

    ## Fitting zero probability model ...

    ## Estimating densities ...

    ## Configuring design ...

    ## Simulating data ...

    ##  ...1 of 1

``` r
sim.data.sc1 <- sim.data.sc[[1]]
class(sim.data.sc1)
```

    ## [1] "SingleCellExperiment"
    ## attr(,"package")
    ## [1] "SingleCellExperiment"

``` r
head(counts(sim.data.sc1)[, seq_len(5)])
```

    ##        Sample_1 Sample_2 Sample_3 Sample_4 Sample_5
    ## Gene_1        0        0       48        0        0
    ## Gene_2       81        7        0        0        0
    ## Gene_3        0       16        0        8      159
    ## Gene_4        0        0        0        0        8
    ## Gene_5       22        0        2        0        0
    ## Gene_6        0        0        0        0       29

``` r
colData(sim.data.sc1)
```

    ## DataFrame with 100 rows and 3 columns
    ##               Batch    Group sim.Lib.Size
    ##            <factor> <factor>    <numeric>
    ## Sample_1          1        1        62456
    ## Sample_2          1        1        64650
    ## Sample_3          1        1       101039
    ## Sample_4          1        1        64104
    ## Sample_5          1        1        60436
    ## ...             ...      ...          ...
    ## Sample_96         1        2        96673
    ## Sample_97         1        2        46853
    ## Sample_98         1        2        58025
    ## Sample_99         1        2        88526
    ## Sample_100        1        2        77400

``` r
rowData(sim.data.sc1)
```

    ## DataFrame with 2000 rows and 2 columns
    ##              DE.ind       source.ID
    ##           <logical>        <factor>
    ## Gene_1         TRUE ENSG00000093072
    ## Gene_2         TRUE ENSG00000169446
    ## Gene_3         TRUE ENSG00000125703
    ## Gene_4         TRUE ENSG00000081320
    ## Gene_5         TRUE ENSG00000106346
    ## ...             ...             ...
    ## Gene_1996     FALSE ENSG00000153767
    ## Gene_1997     FALSE ENSG00000006432
    ## Gene_1998     FALSE ENSG00000168291
    ## Gene_1999     FALSE ENSG00000108774
    ## Gene_2000     FALSE ENSG00000069424
