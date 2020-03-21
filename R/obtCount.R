# Calculates height and mid points of a distribution
# 
# @param Y a vector of gene expression data for a particular gene (in log CPM)
# @param w a numeric value between 0 and 1 or NULL refering the number of classes to be created
# for the outcome data (if NULL the algorithm in graphics::hist() function will be used)
# @param  ... further arguments passed to or from other methods.
# 
# @return a list object contating class mid points, counts, and others
# @examples 
# #
# @importFrom graphics hist
obtCount <- function(Y, w, ...){
  
  if(is.null(w)){
    h   <- hist(Y, plot = FALSE, right = TRUE)
  }
  else if(w>0 & w<1){
    h   <- hist(Y, nclass = round(w*length(Y)), plot = FALSE, right = TRUE)
  }
  s   <- h$breaks
  lls <- s[seq_len(length(s)-1)]
  uls <- s[2:(length(s))]
  S   <- (lls+uls)/2
  Ny <- h$counts
  list(S = S, lls = lls, uls = uls, Ny = Ny, mu.hat = mean(Y), sd.hat = sd(Y)) 
}