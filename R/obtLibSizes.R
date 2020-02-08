# A function to calculate the library sizes in the source data
obtLibSizes <- function(s.data){
  if(is(s.data, "SingleCellExperiment")){ 
    LS <- colSums(counts(s.data))
  }else if(is(s.data, "data.frame") | is(s.data, "matrix")){
    LS <- colSums(s.data)
  }
  return(LS)
}