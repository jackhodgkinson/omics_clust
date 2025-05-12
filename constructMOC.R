# constructMOC.R
constructMOC <- function(data,                    # Input as data frame or list of data frames.
                         ID = NULL,               # ID column for participants
                         parallel = FALSE          # Use parallel processsing. Default is TRUE.
                         )              

  {
    # Load relevant libraries
    library(parallel)

    # Generate MOC for a single dataset
    if (class(data) != "list") {
      data <- as.data.frame(sapply(data, as.factor))
      moc <- do.call(cbind, lapply(data, function(x) model.matrix(~ x - 1)))
      moc <- as.matrix(moc)
      colnames(moc) <- NULL
      names(moc) <- "MOC"  # Ensure proper naming for the single MOC
    }
    
    # Generate MOC for multiple datasets
    else {
      
      # Only keep individuals that are common across all datasets
      if (length(unique(sapply(data, nrow))) != 1) {
        if(!is.null(ID)){
          common_IDs <- Reduce(intersect, lapply(data, function(data) data$ID))
          data <- lapply(data, function(data) data[data$ID %in% common_ids, ])
        }
        else {
          stop("If datasets are not provided with an ID column, then they need
                to have all the same number of rows (and the elements must be in
                the same order in each dataset), otherwise it is impossible
                to match observations from different datasets into the same
                matrix of clusters.")
        }
      }
      
      # With parallel processing
      if (parallel == TRUE) {
      
        # Create empty list for MOC
        moc <- list()
        
        # Create a cluster with available cores
        cl <- makeCluster(detectCores() - 1)
        
        # Export necessary variables and functions
        clusterExport(cl, c("data"))
        
        # Create MOC for each group of data
        moc <- parLapply(cl, 1:length(data), function(i) {
          data[[i]] <- as.data.frame(sapply(data[[i]], as.factor))
          m <- do.call(cbind, lapply(data[[i]], function(x) model.matrix(~ x - 1)))
          m <- as.matrix(m)
          colnames(m) <- NULL
          return(m)
        })
        
        names(moc) <- paste0("MOC - Group", 1:length(moc))
        
        stopCluster(cl)  # Stop the cluster
      }
    
      # Without parallel processing
      else  {
        moc <- list()
        for (i in 1:length(data)) {
          data[[i]] <- as.data.frame(sapply(data[[i]], as.factor))
          moc[[i]] <- do.call(cbind, lapply(data[[i]], function(x) model.matrix(~ x - 1)))
          moc[[i]] <- as.matrix(moc[[i]])
          colnames(moc[[i]]) <- NULL
          names(moc)[[i]] <- paste0("MOC - Group", i)  # Ensure proper naming in the non-parallel block
        }
      }
    }
    return(moc)
  }