#' A function to prepare outputs
#' @param sim.dat The simulated data
#' @param exprmt.design the design 
#' @param DE.ind the differential abundance indicator
#' @param result.format the desired output format
#' @param LL simulated library sizes
#' 
#' @return the data in the desired format
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
prepareSPsimOutputs <- function(sim.dat, exprmt.design, DE.ind, result.format, 
                                LL){
  sim.dat[is.na(sim.dat)] <- 0
  source.ID = rownames(sim.dat)
  rownames(sim.dat) <- paste0("Gene_", seq_len(nrow(sim.dat)))
  colnames(sim.dat) <- paste0("Sample_", seq_len(ncol(sim.dat)))
  col.data <- data.frame("Batch" = exprmt.design$sub.batchs,
             "Group"  = exprmt.design$sub.groups,
             "sim.Lib.Size" = LL,
             row.names = colnames(sim.dat))
  
  row.data <- data.frame("DE.ind" = DE.ind, row.names = rownames(sim.dat),
                         "source.ID" = source.ID)
  sim.data = if(result.format == "SCE"){
     SingleCellExperiment(assays = list(counts = sim.dat),
                          colData = col.data, rowData = row.data)
    } else if(result.format == "list"){
      list("counts" = sim.dat, "colData" = col.data, "rowData" = row.data)
    }
  return(sim.data)
}