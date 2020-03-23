# Check for data validity
checkInputValidity <- function(s.data, group, batch, group.config, batch.config,
                               w, log.CPM.transform, prior.count, 
                              lib.size.params, llStat.thrld){
  # Check for class of missing values in the source data
  if(anyNA(group)){
    stop("The group indicator contains missing values!")
  }
  if(anyNA(batch)){
    stop("The batch indicator contains missing values!")
  }
  if(anyNA(s.data)){
    stop("The source count matrix contains missing values")
  } 
  # Check for group/batch coformation with count matrix
  if(!is.null(group)){
    if(length(group) != ncol(s.data)){
      stop("The length of the group variable is not conformable with the 
      number of columns in the source count matrix!")
    }
  }
  if(!is.null(batch)){
    if(length(group) != ncol(s.data)){
      stop("The length of the batch variable is not conformable with the 
      number of columns in the source count matrix!")
    }
  }
  # Check for the validity of the group/batch configuration
  if(is.null(group) & (length(group.config)!=1)){
   stop("If group is NULL, the length of group.config must equal to 1!")
  }
  if(is.null(batch) & (length(batch.config)!=1)){
    stop("If batch is NULL, the length of batch.config must equal to 1!")
  }
  if(sum(group.config)!=1 | sum(batch.config)!=1){
    stop("The group.cofig and batch.config must sum to 1!")
  }
  # Check for the validity of other inputs
  if(!is.null(group)){
    if(length(group.config) > length(unique(group))){
      stop("The number of groups to be simulated is larger than the number of 
      groups available in the source data!")
    } else if((length(group.config) < length(unique(group))) & length(group.config) !=1){
     stop("The number of groups to be simulated is less than the number of 
      groups available in the source data. Maybe subset your source data OR consider 
      adding a 0 in the group.config vector for the group that you do not want to simulate 
      from. Example group.config=c(0.5, 0.5, 0)!")
    } else if((length(unique(group))!=1) & length(group.config) ==1){
      stop("The source data has", length(unique(group)), "groups but you are
                      simulating only for one group (length(group.config)=1). Hence, data 
                      only from one of the groups will be used!")
    } 
  }else{
    if(length(group.config) > 1){
      stop("You are attempting to simulate",  length(group.config), "groups but 
                   there is no groupping variable (group=NULL) for the source data. \n")
    }
  } 
  if(!is.null(batch)){
    if(length(batch.config) != length(unique(batch))){
      stop("The number of batchs to be simulated is not equal to the number of batches 
      available in the source data!")
    }
  }else{
    if(length(batch.config) > 1){
      stop("You are attempting to simulate",  length(batch.config), "batches but 
                     there is no batch indicator (batch=NULL) for the source data. \n")
    }  
  }   
  if(!is.null(lib.size.params) & (length(lib.size.params) != 2 | is.null(names(lib.size.params)))){
    stop("The log-normal parameters for the distribution of library sizes must be submitted in a named vector of size 2. 
             Example, lib.size.params = c(meanlog=10, sdlog=0.2). See also ?rlnorm()")
  }
  if(!is.null(w) && (w < 0 || w > 1)){
    stop("w should be NULL or any value between 0 and 1 excluding 0 and 1")
  }
  if(prior.count < 0){
    stop("Prior count!")
  }
  if(llStat.thrld <0){
    stop("Likelihood ratio test statistic threshold should be non-negative")
  }
  invisible()
}

