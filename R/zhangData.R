#' Bulk RNA-seq data from Zhang et al study.
#'
#' 498 neuroblastoma tumors, was retrieved from Zhang et al. (GEO accession number GSE49711). 
#' In short, unstranded poly(A)+ RNA sequencing was performed on the HiSeq 2000 
#' instrument (Illumina). Paired-end reads with a length of 100 nucleotides 
#' were obtained. To quantify the full transcriptome, raw fastq files 
#' were processed with Kallisto v0.42.4 (index build with GRCh38-Ensembl v85). 
#' The pseudo-alignment tool Kallisto was chosen above other quantification methods 
#' as it is performing equally good but faster. For this study, a subset of 
#' 172 patients with high-risk disease were selected, forming two groups: the MYCN amplified (n1 = 91) 
#' and MYCN non-amplified (n2 = 81) tumors.
#'
#' @format A list object
#' @references Zhang W, Yu Y, Hertwig F, Thierry-Mieg J, Zhang W, Thierry-Mieg D, Wang J, Furlanello C, Devanarayan V, Cheng J, et al. Comparison of RNA-seq and microarray-based models for clinical endpoint prediction. Genome Biol. 2015;16(133) https://doi.org/10.1186/s13059-015-0694-1
#' \describe{
#'   \item{counts}{gene counts}
#'   \item{group}{MYCN (0 for MYCN non-amplified and 1 for MYCN amplified)}
#' }
#' @source \url{https://doi.org/10.1186/s13059-015-0694-1}
"zhang.data"