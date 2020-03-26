#' Generate library sizes from log-normal
#'
#' @param fit.ln 
#' @param exprmt.design 
#'
#' @return The generated libray sizes per batch and group
genLibSizes = function(fit.ln, exprmt.design){
  Seq = seq_along(n.group)
  n.batch = rowSums(exprmt.design$exprmt.config)
  lapply(unique(exprmt.design$sub.batchs), function(b){
    rlnorm(n.batch[b], fit.ln[[b]][["meanlog"]], fit.ln[[b]][["sdlog"]])
  })
}