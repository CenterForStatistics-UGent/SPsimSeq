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

