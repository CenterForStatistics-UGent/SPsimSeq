#' Generate library sizes from log-normal
#'
#' @param n.batch 
#' @param n.group 
#' @param fit.ln 
#' @param config.mat 
#'
#' @return The generated libray sizes per batch and group
genLibSizes = function(n.batch, n.group, fit.ln, config.mat){
  Seq = seq_along(n.group)
  lapply(seq_along(fit.ln), function(b){
    L.b.pred <- rlnorm(n.batch[b], fit.ln[[b]][["meanlog"]], fit.ln[[b]][["sdlog"]])
    ## Split into the groups
    gr <- rep(Seq, config.mat[b, ])
    split(L.b.pred, gr) 
  })
}