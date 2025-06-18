# constructMOC.R
# Load numCores.R
source("numCores.R")

# Function
constructMOC <- function(classification,           # Input as data frame or list of data frames.
                         ID = NULL,                # ID column for participants
                         parallel_process = FALSE  # Use parallel processsing. Default is TRUE.
                         )              

  {
    # Load relevant libraries
    library(parallel)
  
    # Check if input is a matrix when not a list
    if (!is.list(classification) && is.matrix(classification)) {
      stop("Input 'classification' must be a data frame, not a matrix. Use as.data.frame() before passing it.")
    }

    # Generate MOC for a single dataset
    if (class(classification) != "list") {
      classification <- as.data.frame(sapply(classification, as.factor))
      moc <- do.call(cbind, lapply(classification, function(x) model.matrix(~ x - 1)))
      moc <- as.matrix(moc)
      colnames(moc) <- NULL
      names(moc) <- "MOC"  # Ensure proper naming for the single MOC
    }
    
    # Generate MOC for multiple datasets
    else {
      
      # Only keep individuals that are common across all datasets
      if (length(unique(sapply(classification, nrow))) != 1) {
        if(!is.null(ID)){
          common_IDs <- Reduce(intersect, lapply(classification, function(classification) classification$ID))
          classification <- lapply(classification, function(classification) classification[classification$ID %in% common_ids, ])
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
      if (parallel_process) {
      
        # Create empty list for MOC
        moc <- list()
        
        # If using Windows
        if (.Platform$OS.type == "windows"){
        
          # Create a cluster with available cores
          cl <- makeCluster(n_cores)
          
          # Export necessary variables and functions
          clusterExport(cl, c("data"))
          
          # Create MOC for each group of data
          moc <- parLapply(cl, 1:length(classification), function(i) {
            classification[[i]] <- as.data.frame(sapply(classification[[i]], as.factor))
            m <- do.call(cbind, lapply(classification[[i]], function(x) model.matrix(~ x - 1)))
            m <- as.matrix(m)
            colnames(m) <- NULL
            return(m)
          })
          
          # End the cluster
          stopCluster(cl)
        
        } else {
          moc <- mclapply(1:length(classification), function(i) {
            data[[i]] <- as.data.frame(sapply(classification[[i]], as.factor))
            m <- do.call(cbind, lapply(classification[[i]], function(x) model.matrix(~ x - 1)))
            m <- as.matrix(m)
            colnames(m) <- NULL
            return(m)
          }, mc.cores = n_cores)
        } 
        names(moc) <- paste0("MOC - Group", 1:length(moc))
      }
    
      # Without parallel processing
      else  {
        moc <- list()
        for (i in 1:length(classification)) {
          classification[[i]] <- as.data.frame(sapply(classification[[i]], as.factor))
          moc[[i]] <- do.call(cbind, lapply(classification[[i]], function(x) model.matrix(~ x - 1)))
          moc[[i]] <- as.matrix(moc[[i]])
          colnames(moc[[i]]) <- NULL
          names(moc)[[i]] <- paste0("MOC - Group", i)  
        }
      }
    }
    return(moc)
  }