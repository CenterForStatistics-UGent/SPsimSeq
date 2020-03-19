#' Fast fit Poisson regression
#'
#' @param Ny vector of counts
#' @param x regressor
#' @param degree degree of the polynomial
#' @param offset offset
#' @importFrom stats glm.fit poisson
#' @return see glm.fit
fitPoisGlm = function(Ny, x, degree, offset){
  desMat = buildXmat(x, degree+1)
  glm.fit(y = Ny, x = desMat, family = poisson(), offset = offset)
}