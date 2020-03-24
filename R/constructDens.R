#' Construct the cumulative density
#'
#' @param densList.ii the estimated density parameters
#' @param exprmt.config experiment configuration
#' @param DE.ind.ii a boolean, is the gene to be DE?
#'
#' @return The cumulative density
constructDens = function(densList.ii, exprmt.config, DE.ind.ii){
    lapply(seq_along(exprmt.config$sub.batchs), function(i){
      batch = exprmt.config$sub.batchs[[i]]
      dl = if(DE.ind.ii){
        densList.ii[[batch]][[exprmt.config$sub.group[[i]]]]
      } else {
        densList.ii[[batch]]
      }
      Coef = dl$coef
      s.mat <- buildXmat(dl$mids, nc = length(Coef))
      gy = dl$g0*dl$n*exp(s.mat %*% Coef)
      id = is.infinite(gy)
      if(any(id)){
        gy[id] = max(gy[!id])
      }
      Gy = cumsum(gy)/sum(gy)
      return(Gy)
    })
  }
