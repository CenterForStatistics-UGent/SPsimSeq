# A function to prepare outputs
prepareSPsimOutputs <- function(sim.dat, n.batch, n.group, config.mat, LL, DE.ind, sel.genes,
                                result.format, log.CPM.transform){
  sim.count <- sim.dat
  sim.count[is.na(sim.count)] <- 0
  rownames(sim.count) <- paste0("Gene_", seq_len(nrow(sim.count)))
  colnames(sim.count) <- paste0("Sample_", seq_len(ncol(sim.count)))
  
  # col.data <- data.frame("Batch"=rep(rep(seq_len(length(n.batch)), times=length(n.group)), config.mat),
  #                        Group=rep(rep(seq_len(length(n.group)), each=length(n.batch)), config.mat), 
  #                        sim.Lib.Size=do.call("c", do.call("c", LL)), row.names = colnames(sim.count))
  if(log.CPM.transform){
    sim.Lib.Size = do.call("c", do.call("c", LL))
  }else{
    sim.Lib.Size = NA
  }
  col.data <- data.frame("Batch" = do.call("c",lapply(seq_along(n.batch), 
                                          function(i) rep(i, n.batch[i]))),
             "Group"  = do.call("c",lapply(seq_along(n.group), 
                                           function(i) rep(i, n.group[i]))),
             "sim.Lib.Size" = sim.Lib.Size,
             row.names = colnames(sim.count))
  
  row.data <- data.frame("DE.ind"=DE.ind, "source.ID"=sel.genes, "row.names" = rownames(sim.count))
  
  # SPsim.est.densities <- lapply(seq_len(length(sim.list)), function(ii){ 
  #   list(SPsim.dens.hat=sim.list[[ii]]$SPsim.dens.hat,  
  #        norm.carrier.dens = sim.list[[ii]]$norm.carrier.dens,
  #        SPsim.dens.parms=sim.list[[ii]]$SPsim.dens.parms, 
  #        zero.prob.params=sim.list[[ii]]$zero.prob.params,
  #        source_gene.ID = sel.genes[ii])
  # })
  # names(SPsim.est.densities) <- rownames(sim.count)
  
  
  if(result.format == "SCE"){
    sim.data <- SingleCellExperiment(assays=list(counts=sim.count),
                                     colData=col.data, rowData=row.data)
    #sim.data <- sim.data[rowSums(counts(sim.data)>0)>0,]
    return(sim.data)
  }
  else if(result.format == "list"){
    #keep <- rowSums(sim.count>0)>0
    sim.data <- list("counts"=sim.count, "colData"=col.data, "rowData"=row.data)
    # sim.data <- list("counts"=sim.count[keep,], "colData"=col.data, "rowData"=row.data[keep,],
    #                  SPsim.est.densities=SPsim.est.densities[keep])
    return(sim.data)
  }
}