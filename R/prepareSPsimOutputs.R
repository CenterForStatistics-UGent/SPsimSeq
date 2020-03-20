# A function to prepare outputs
prepareSPsimOutputs <- function(sim.dat, n.batch, n.group, LL, DE.ind, sel.genes,
                                result.format, log.CPM.transform){
  sim.count <- t(sapply(sim.dat, function(x) x[,"y.new"]))
  sim.count[is.na(sim.count)] <- 0
  rownames(sim.count) <- paste0("Gene_", seq_len(nrow(sim.count)))
  colnames(sim.count) <- paste0("Sample_", seq_len(ncol(sim.count)))
  if(log.CPM.transform){
    sim.Lib.Size = do.call("c", do.call("c", LL))
  }else{
    sim.Lib.Size = NA
  }
  col.data <- data.frame("Batch" = sim.dat[[1]]$Batch,
             "Group"  = sim.dat[[1]]$Group,
             "sim.Lib.Size" = sim.Lib.Size,
             row.names = colnames(sim.count))
  
  row.data <- data.frame("DE.ind"=DE.ind, "source.ID"=sel.genes, "row.names" = rownames(sim.count))
  if(result.format == "SCE"){
    sim.data <- SingleCellExperiment(assays=list(counts=sim.count),
                                     colData=col.data, rowData=row.data)
    return(sim.data)
  }
  else if(result.format == "list"){
    sim.data <- list("counts"=sim.count, "colData"=col.data, "rowData"=row.data)
    return(sim.data)
  }
}