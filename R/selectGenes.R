#' Sample genes from candidate genes
#'
#' @param pDE fraction of genes to be made DE
#' @param exprmt.design the experiment design
#' @param n.genes the total number of genes required
#' @param null.genes0,nonnull.genes0 Candidate null and non-null genes
#'
#' @return a vector of selected genes
selectGenes <- function(pDE, exprmt.design, n.genes, null.genes0, nonnull.genes0){
  if(pDE>0 && length(unique(exprmt.design$sub.groups))>1){
    null.genes = sample(null.genes0, (1-pDE)*(n.genes), 
                        replace = (1-pDE)*n.genes > length(null.genes0))
    nonnull.genes = if(length(nonnull.genes0)>0){
          sample(nonnull.genes0, pDE*n.genes, 
                 replace = pDE*n.genes > length(nonnull.genes0))
    } else {NULL}
    sel.genes <- c(nonnull.genes, null.genes)
  } else{
    sel.genes <- sample(null.genes0, n.genes, 
                        replace = n.genes <= length(null.genes0))
  } 
  names(sel.genes) = sel.genes
  sel.genes
}