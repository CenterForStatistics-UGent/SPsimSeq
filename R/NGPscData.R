#' Neuroblastoma NGP cells single-cell RNA-seq.
#'
#'  It was retrieved from [1] (GEO accession GSE119984): This dataset is 
#'  generated for a cellular perturbation experiment on the C1 instrument 
#'  (SMARTer protocol) [1]. This total RNA-seq dataset contains 
#'  83 NGP neuroblastoma cells, of which 31 were treated with 
#'  nutlin-3 and the other 52 cells were treated with vehicle (controls). 
#' 
#' @docType data
#' 
#' @format A SingleCellExperiment object
#' 
#' @keywords dataset
#' 
#' @references 1 - Verboom, K., Everaert, C., Bolduc, N., Livak, K. J., Yigit, N., Rombaut, D., ... & Speleman, F. (2019). SMARTer single cell total RNA sequencing. Nucleic Acids Research, 47(16), e93-e93.
#' 
#' \describe{
#'   \item{SingleCellExperiment}{counts + gene info + cell infro} 
#' }
#'  
#' @source GEO accession GSE119984
#' @examples 
#' data("scNGP.data")
#' scNGP.data
"scNGP.data"

