# A function to prepare outputs
prepareSPsimOutputs <- function(sim.list, n.batch, n.group, config.mat, LL, DE.ind, sel.genes,
                                result.format){
  sim.count <- do.call(rbind, lapply(sim.list, function(x){x$Y.star}))
  sim.count[is.na(sim.count)] <- 0
  rownames(sim.count) <- paste0("Gene_", seq_len(nrow(sim.count)))
  colnames(sim.count) <- paste0("Sample_", seq_len(ncol(sim.count)))
  
  col.data <- data.frame(Batch=rep(rep(seq_len(length(n.batch)), times=length(n.group)), config.mat),
                         Group=rep(rep(seq_len(length(n.group)), each=length(n.batch)), config.mat), 
                         sim.Lib.Size=do.call("c", do.call("c", LL)), row.names = colnames(sim.count))
  row.data <- data.frame(DE.ind=DE.ind, source.ID=sel.genes, row.names = rownames(sim.count))
  
  SPsim.est.densities <- lapply(sim.list, function(x){
    list(SPsim.dens.hat=x$SPsim.dens.hat,  norm.carrier.dens = x$norm.carrier.dens)
  })
  
  if(result.format == "SCE"){
    sim.data <- SingleCellExperiment(assays=list(counts=sim.count),
                                     colData=col.data, rowData=row.data,
                                     metadata=list(SPsim.est.densities=SPsim.est.densities))
    sim.data 
  }
  else if(result.format == "list"){
    sim.data <- list(counts=sim.count, colData=col.data, rowData=row.data,
                     SPsim.est.densities=SPsim.est.densities)
    sim.data
  }
}