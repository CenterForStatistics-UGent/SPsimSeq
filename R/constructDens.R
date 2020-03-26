#' Construct the cumulative density
#'
#' @param densList.ii the estimated density parameters
#' @param exprmt.design experiment configuration
#' @param DE.ind.ii a boolean, is the gene to be DE?
#'
#' @return The cumulative density
constructDens = function(densList.ii, exprmt.design, DE.ind.ii){
    lapply(seq_along(exprmt.design$sub.batchs), function(i){
      batch = exprmt.design$sub.batchs[[i]]
      dl = if(DE.ind.ii){
        densList.ii[[batch]][[exprmt.design$sub.group[[i]]]]
      } else {
        densList.ii[[batch]]
      }
      Coef = dl$coef
      if(is.null(Coef)){
        gy = dl$g0*dl$n
      } else {
        s.mat <- buildXmat(dl$mids, nc = length(Coef))
        gy = dl$g0*dl$n*exp(s.mat %*% Coef)
        id = is.infinite(gy)
        if(any(id)){
          gy[id] = max(gy[!id])
        }
      }
      Gy = cumsum(gy)/sum(gy)
      return(list(Gy = Gy, breaks = dl$breaks))
    })
  }
