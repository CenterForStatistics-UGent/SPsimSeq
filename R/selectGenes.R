# selectes genes among candidate genes
selectGenesSim <- function(pDE, group, n.genes, null.genes0, nonnull.genes0,
                           group.config, ...){
  if(pDE>0 & !is.null(group) & length(group.config)>1){
    if((1-pDE)*n.genes <= length(null.genes0)){
      #set.seed(setSimContol$seed.sample.genes[h])
      null.genes     <- sample(null.genes0, (1-pDE)*(n.genes), replace = FALSE)
    }
    else{
      message("Note: The number of null genes (not DE) in the source data is ", length(null.genes0), 
              " and the number of null genes required to be included in the simulated data is ",
              round((1-pDE)*(n.genes)), ". Therefore, candidiate null genes are sampled with replacement.")
      #set.seed(setSimContol$seed.sample.genes[h])
      null.genes     <- sample(null.genes0, (1-pDE)*(n.genes), replace = TRUE)
    }
    
    
    if(pDE*n.genes <= length(nonnull.genes0)){
      #set.seed(setSimContol$seed.sample.genes[h]+1)
      nonnull.genes  <- sample(nonnull.genes0, pDE*n.genes, replace = FALSE)
    }
    else{
      message("Note: The number of DE genes detected in the source data is ", length(nonnull.genes0), 
              " and the number of DE genes required to be included in the simulated data is ",
              round(pDE*n.genes), ". Therefore, candidiate DE genes are sampled with replacement.")
      #set.seed(setSimContol$seed.sample.genes[h]+1)
      nonnull.genes  <- sample(nonnull.genes0, pDE*n.genes, replace = TRUE)
    }
    #set.seed(NULL)
    
    sel.genes <- c(nonnull.genes, null.genes)
    DE.ind <- ifelse(sel.genes %in% null.genes, 0, 1)
    names(DE.ind) <- sel.genes
    
  }
  else{
    if(n.genes <= length(null.genes0)){
      #set.seed(setSimContol$seed.sample.genes[h])
      sel.genes     <- sample(null.genes0, n.genes, replace = FALSE)
    }
    else{
      #set.seed(setSimContol$seed.sample.genes[h])
      sel.genes     <- sample(null.genes0, n.genes, replace = TRUE)
    } 
    DE.ind <- rep(0, length(sel.genes))
    names(DE.ind) <- sel.genes
  } 
  
  list(DE.ind=DE.ind, sel.genes=sel.genes)
}