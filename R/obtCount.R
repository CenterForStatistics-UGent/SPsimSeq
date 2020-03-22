#' Calculates height and mid points of a distribution
#' 
#' @param Y a vector of gene expression data for a particular gene (in log CPM)
# 
#' @return a list object contating class mid points, counts, and others
#' @importFrom graphics hist
obtCount <- function(Y){
  h   <- hist(Y, plot = FALSE, right = TRUE)
  s   <- h$breaks
  lls <- s[seq_len(length(s)-1)]
  uls <- s[2:(length(s))]
  S   <- (lls+uls)/2
  Ny <- h$counts
  list(S = S, lls = lls, uls = uls, Ny = Ny, mu.hat = mean(Y), sig.hat = sd(Y)) 
}