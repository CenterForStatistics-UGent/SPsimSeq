#' Generate library sizes from log-normal
#'
#' @param fit.ln 
#' @param exprmt.design 
#'
#' @return The generated libray sizes per batch and group
#' @importFrom stats rlnorm
genLibSizes = function(fit.ln, exprmt.design){
  n.batch = rowSums(exprmt.design$exprmt.config)
  lapply(unique(exprmt.design$sub.batchs), function(b){
    rlnorm(n.batch[b], fit.ln[[b]][["meanlog"]], fit.ln[[b]][["sdlog"]])
  })
}