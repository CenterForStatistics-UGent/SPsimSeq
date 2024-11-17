#' Calculates height and mid points of a distribution
#' 
#' @param Y a vector of gene expression data for a particular gene (in log CPM)
#' @param w a numeric value between 0 and 1 or NULL refering the number of classes to be created
# for the outcome data (if NULL the algorithm in graphics::hist() function will be used)
# 
#' @return a list object contating class breaks, mid points and counts
#' @importFrom graphics hist
obtCount <- function(Y, w){
  nclass=if(is.null(w)) NULL else round(w*length(Y))
  h = hist(Y, plot = FALSE, nclass=nclass) 
  list(breaks = h$breaks, mids = h$mids, counts = h$counts) 
}
