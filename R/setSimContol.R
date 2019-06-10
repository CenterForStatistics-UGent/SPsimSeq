# sets seed for reproducible siulations
setSimContol <- function(seed, n.sim, n.genes){
  
  set.seed(seed)
  seed.batch.pars <- sample(1:1e9, n.genes)
  
  set.seed(seed)
  seed.sample.genes <- sample(1:1e9, n.sim)
  
  set.seed(seed)
  seed.sample.Y   <- sample(1:1e9, n.genes)
  
  set.seed(seed)
  seed.sample.single <- sample(1:1e4, 1) 
  
  
  list(seed.batch.pars=seed.batch.pars, seed.sample.Y=seed.sample.Y, 
       seed.sample.genes=seed.sample.genes, seed.sample.single=seed.sample.single)
} 
