# Check for data validity
checkInputValidity <- function(s.data, group, batch, group.config, batch.config){
  # Check for class of missing values in the source data
  if(anyNA(group)){
    val2 <- -1
    mes2 <-  "    The group indicator missing value. \n"
  }else{val2 <- 0; mes2 <- "none"}
  if(any(is.na(batch))){
    val3  <- -1
    mes3 <-  "The batch indicator missing value. \n"
  }else{val3 <- 0; mes3 <- "none"}
  
  if(is(s.data, "SingleCellExperiment")){
    if(any(is.na(counts(s.data)))){
      val4  <- -1
      mes4 <- "The source count matrix contains missing value. \n"
    }else{val4 <- 0; mes4 <- "none"}
    
  }else{
    if(any(is.na(s.data))){
      val4  <- -1
      mes4 <-  "    The source count matrix contains missing value. \n"
    }else{val4 <- 0; mes4 <- "none"}  
  }
  
  # Check for group/batch coformation with count matrix
  if(!is.null(group)){
    if(length(group) != ncol(s.data)){
      val5  <- -1
      mes5 <-  "The length of the group variable is not conformable with the 
      number of columns in the source count matrix. \n"
    }else{val5 <- 0; mes5 <- "none"} 
  }else{val5 <- 0; mes5 <- "none"} 
  
  if(!is.null(batch)){
    if(length(batch) != ncol(s.data)){
      val6  <- -1
      mes6 <-  "The length of the batch variable is not conformable with the 
      number of columns in the source count matrix. \n"
    }else{val6 <- 0; mes6 <- "none"} 
  }else{val6 <- 0; mes6 <- "none"} 
  # Check for the validity of the group/batch configuration
  if(is.null(group) & (length(group.config)!=1)){
    val9  <- -1
    mes9 <-  "    If group is NULL, the length of group.config must equal to 1. \n"
  }else{val9 <- 0; mes9 <- "none"}
  
  if(is.null(batch) & (length(batch.config)!=1)){
    val10  <- -1
    mes10 <-  "    If batch is NULL, the length of batch.config must equal to 1. \n"
  }else{val10 <- 0; mes10 <- "none"}
  
  if(sum(group.config)!=1 | sum(batch.config)!=1){
    val11  <- -1
    mes11 <-  "    The group.cofig and batch.config must sum to 1. \n"
  }else{val11 <- 0; mes11 <- "none"}
  
  
  
  # Check for the validity of other inputs
  if(!is.null(group)){
    if(length(group.config) > length(unique(group))){
      val12  <- -1
      mes12 <-  "The number of groups to be simulated is larger than the number of 
      groups available in the source data. \n"
    }
    
    else if((length(group.config) < length(unique(group))) & length(group.config) !=1){
      val12  <- -1
      mes12 <-  "The number of groups to be simulated is less than the number of 
      groups available in the source data. Maybe subset your source data OR consider 
      adding a 0 in the group.config vector for the group that you do not want to simulate 
      from. Example group.config=c(0.5, 0.5, 0). \n"
    } 
    else if((length(unique(group))!=1) & length(group.config) ==1){
      val12  <- 1 
      mes12 <-  paste("The source data has", length(unique(group)), "groups but you are
                      simulating only for one group (length(group.config)=1). Hence, data 
                      only from one of the groups will be used. \n")
    } 
    else{val12 <- 0; mes12 <- "none"} 
  }else{
    if(length(group.config) > 1){
      val12  <- -1
      mes12 <- paste("You are attempting to simulate",  length(group.config), "groups but 
                   there is no groupping variable (group=NULL) for the source data. \n")
    } 
    else{val12 <- 0; mes12 <- "none"}
  } 
  
  if(!is.null(batch)){
    if(length(batch.config) != length(unique(batch))){
      val13  <- 1 
      mes13 <- "The number of batchs to be simulated is not equal to the number of batches 
      available in the source data. \n"
    }else{val13 <- 0; mes13 <- "none"}  
  }else{
    if(length(batch.config) > 1){
      val13  <- -1
      mes13 <- paste("You are attempting to simulate",  length(batch.config), "batches but 
                     there is no batch indicator (batch=NULL) for the source data. \n")
    }else{val13 <- 0; mes13 <- "none"}   
  }   
  
  val <- sapply(paste0("val", seq_len(13)[-c(7, 8)]), function(x) eval(parse(text = x)))
  mes <- sapply(paste0("mes", seq_len(13)[-c(7, 8)]), function(x) eval(parse(text = x)))
  
  list(val=val, mes=mes)
}

