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
BiocManager::install("SPsimSeq")
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

    ## Selecting candidate DE genes ...

    ## Note: The number of DE genes detected in the source data is 66 and the number of DE genes required to be included in the simulated data is 200. Therefore, candidiate DE genes are sampled with replacement.

    ## Estimating densities ...

    ## Configuring design ...

    ## Simulating data ...

    ##  ...1 of 1

``` r
sim.data.bulk1 <- sim.data.bulk[[1]]                              
head(sim.data.bulk1$counts[, seq_len(5)])  # count data
```

    ##        Sample_1 Sample_2 Sample_3 Sample_4 Sample_5
    ## Gene_1      124        1       57      134      108
    ## Gene_2        1        3       44        0       91
    ## Gene_3        0        3       75        1      126
    ## Gene_4       11       10       24       10        6
    ## Gene_5        0        1        0        0        2
    ## Gene_6      187       31      332       27       73

``` r
head(sim.data.bulk1$colData)               # sample info
```

    ##          Batch Group sim.Lib.Size
    ## Sample_1     1     1      1499454
    ## Sample_2     1     1      1403642
    ## Sample_3     1     1      1395046
    ## Sample_4     1     1      1047452
    ## Sample_5     1     1      1428862
    ## Sample_6     1     1      1365620

``` r
head(sim.data.bulk1$rowData)               # gene info
```

    ##        DE.ind      source.ID
    ## Gene_1   TRUE RP11-1006G14.4
    ## Gene_2   TRUE  RP11-706O15.7
    ## Gene_3   TRUE  RP11-706O15.7
    ## Gene_4   TRUE       USP2-AS1
    ## Gene_5   TRUE          XAGE5
    ## Gene_6   TRUE          CAPS2

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

    ## Note: The number of DE genes detected in the source data is 42 and the number of DE genes required to be included in the simulated data is 200. Therefore, candidiate DE genes are sampled with replacement.

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
    ## Gene_1       21        7        0        2        1
    ## Gene_2       19        0        1        0        0
    ## Gene_3       24        0        2       50        0
    ## Gene_4        0       60        0       12       35
    ## Gene_5        0        7        0        8       16
    ## Gene_6        5       16       11        8       41

``` r
colData(sim.data.sc1)
```

    ## DataFrame with 100 rows and 3 columns
    ##               Batch    Group sim.Lib.Size
    ##            <factor> <factor>    <numeric>
    ## Sample_1          1        1        69309
    ## Sample_2          1        1        37024
    ## Sample_3          1        1        81538
    ## Sample_4          1        1       134070
    ## Sample_5          1        1        67218
    ## ...             ...      ...          ...
    ## Sample_96         1        2       122978
    ## Sample_97         1        2        86100
    ## Sample_98         1        2       107930
    ## Sample_99         1        2        59519
    ## Sample_100        1        2        67167

``` r
rowData(sim.data.sc1)
```

    ## DataFrame with 2000 rows and 2 columns
    ##              DE.ind       source.ID
    ##           <logical>        <factor>
    ## Gene_1         TRUE ENSG00000174482
    ## Gene_2         TRUE ENSG00000087586
    ## Gene_3         TRUE ENSG00000164306
    ## Gene_4         TRUE ENSG00000145386
    ## Gene_5         TRUE ENSG00000105613
    ## ...             ...             ...
    ## Gene_1996     FALSE ENSG00000274267
    ## Gene_1997     FALSE ENSG00000143126
    ## Gene_1998     FALSE ENSG00000175305
    ## Gene_1999     FALSE ENSG00000109099
    ## Gene_2000     FALSE ENSG00000237017
