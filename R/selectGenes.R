# selectes genes among candidate genes
selectGenes <- function(pDE, group, n.genes, null.genes0, nonnull.genes0,
                           group.config, ...){
  if(pDE>0 & !is.null(group) & length(group.config)>1){
    null.genes = if((1-pDE)*n.genes <= length(null.genes0)){
      sample(null.genes0, (1-pDE)*(n.genes), replace = FALSE)
    } else{
      sample(null.genes0, (1-pDE)*(n.genes), replace = TRUE)
    }
    
    nonnull.genes = if(length(nonnull.genes0)>0){
      if(pDE*n.genes <= length(nonnull.genes0)){
          sample(nonnull.genes0, pDE*n.genes, replace = FALSE)
      }
      else {
          sample(nonnull.genes0, pDE*n.genes, replace = TRUE)
      }
    } else {NULL}
    sel.genes <- c(nonnull.genes, null.genes)
  } else{
      sel.genes <- sample(null.genes0, n.genes, 
                          replace = n.genes <= length(null.genes0))
  } 
  names(sel.genes) = sel.genes
  sel.genes
}