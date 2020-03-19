# Configure experiment
# 
# @param batch.config a numerical vector for the marginal fraction of samples in each batch. 
# The number of batches to be simulated is equal to the size of the vector.
# All values must sum to 1.
# @param group.config a numerical vector for the marginal fraction of samples in each group.
# The number of groups to be simulated is equal to the size of the vector. All values must sum to 1.
# @param tot.samples total number of samples to be simulated.
# 
# @return a list object contating the number of groups and batches to be simukated, 
# and the experiment configurartion 
# @examples 
# \donttest{
# #---- a design with a total of 10 samples/cells from 1 batch and 1 group
# expriment.config(batch.config=1, group.config=1, tot.samples=10)
# 
# #---- a design with a total of 20 samples/cells from 1 group and 2 batchs with 
# # batch 1 has 15 samples/cells and batch 2 has 5
# expriment.config(batch.config=c(15/20, 5/20), group.config=1, tot.samples=20)
# 
# #---- a design with a total of 20 samples/cells from 1 batch and 2 groups with 
# # group 1 has 10 samples/cells and batch 2 has 10
# expriment.config(batch.config=1, group.config=c(0.5, 0.5), tot.samples=20)
# 
# #---- a design with a total of 30 samples/cells from 2 groups with group 1 has 15 samples 
# # and group 2 has 15, and  three batchs with batch 1,2, and 3 have 5, 10, and 15 samples/cells, 
# # respectively.
# expriment.config(batch.config=c(5/30, 10/30, 15/30), group.config=c(0.5, 0.5), tot.samples=30)
# }
configExperiment <- function(batch.config, group.config, tot.samples = 150){
  if(is.null(names(batch.config))){
    names(batch.config) <- paste0("Batch_",seq_along(batch.config))
  }
  if(is.null(names(group.config))){
    names(group.config) <- paste0("Group_",seq_along(group.config))
  }
  
  exprmt.config <- round(tcrossprod(batch.config, group.config)*tot.samples)
  rownames(exprmt.config) <- names(batch.config)
  colnames(exprmt.config) <- names(group.config)
  
  n.batch <- rowSums(exprmt.config)
  n.group <- colSums(exprmt.config)
  
  list(exprmt.config=exprmt.config, n.batch=n.batch, n.group=n.group)
}