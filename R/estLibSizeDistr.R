#' Estimate log-normal distribution for the library sizes 
#'
#' @param LS observed library sizes
#' @param batch batches
#' 
#' @return Estimated log-normal parameter library sizes
#' @importFrom fitdistrplus fitdistr
estLibSizeDistr = function(LS, batch){
  tapply(LS, batch, function(L.b){
      fitdist(L.b, distr = "lnorm")$estimate 
  })
}