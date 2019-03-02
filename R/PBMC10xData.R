#' Single cell RNA-seq data from 10x genomics
#'
#' Contains 2700 single cells sequenced on an Illumina NextSeq 500 using UMIs. The data is
#' generated using the 10x Genomics Chromium V2 protocol.
#'
#' @format A SingleCellExperiment object
#' \describe{
#'   \item{SingleCellExpeiment}{contains the count data, and gene and cell level 
#'   annotations are stored in rowData and colData slots, respectively.} 
#' }
#' @source \url{https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz}
"PBMC.10x.data"