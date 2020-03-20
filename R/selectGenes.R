# selectes genes among candidate genes
selectGenes <- function(pDE, group, n.genes, null.genes0, nonnull.genes0,
                           group.config, ...){
  if(pDE>0 & !is.null(group) & length(group.config)>1){
    if((1-pDE)*n.genes <= length(null.genes0)){
      null.genes     <- sample(null.genes0, (1-pDE)*(n.genes), replace = FALSE)
    }
    else{
      null.genes     <- sample(null.genes0, (1-pDE)*(n.genes), replace = TRUE)
    }
    
    if(length(nonnull.genes0)>0){
      if(pDE*n.genes <= length(nonnull.genes0)){
        nonnull.genes  <- sample(nonnull.genes0, pDE*n.genes, replace = FALSE)
      }
      else {
        nonnull.genes  <- sample(nonnull.genes0, pDE*n.genes, replace = TRUE)
      }
    }
    else{ 
      nonnull.genes <- NULL
    }
    sel.genes <- c(nonnull.genes, null.genes)
  }  else{
      sel.genes <- sample(null.genes0, n.genes, 
                          replace = n.genes <= length(null.genes0))
  } 
  names(sel.genes) = sel.genes
  sel.genes
}