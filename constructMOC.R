# constructMOC.R
# Load numCores.R
source("numCores.R")

# Function
constructMOC <- function(data,                    # Input as data frame or list of data frames.
                         ID = NULL,               # ID column for participants
                         parallel_process = FALSE         # Use parallel processsing. Default is TRUE.
                         )              

  {
    # Load relevant libraries
    library(parallel)
  
    # # Construct empty similarity matrix 
    # sim_matrix <- matrix(0, nrow = ncol(classification), ncol = ncol(classification),
    #                    dimnames = list(colnames(classification), colnames(classification)))
    # 
    # # Set seed 
    # set.seed(seed)
    # 
    # if (parallel) {
    #   
    #   # Obtain number of cores 
    #   n_cores <- numCores()
    #   
    #   if (.Platform$OS.type == "windows") {
    #     cl2 <- makeCluster(n_cores)
    #     clusterExport(cl2, ls(envir = environment()), envir = environment())
    #     
    #     # Parallelised computation of similarity matrix
    #     sim_mat <- parLapply(cl2, 1:ncol(classification), function(i) {
    #       sapply(1:ncol(classification), function(j) {
    #         adjustedRandIndex(classification[[i]], classification[[j]])
    #       })
    #     })
    #     
    #     # Stop the parallel cluster after the work is done
    #     stopCluster(cl2)
    #     
    #   } else {
    #     sim_mat <- mclapply(1:ncol(classification), function(i) {
    #       sapply(1:ncol(classification), function(j) {
    #         adjustedRandIndex(classification[[i]], classification[[j]])
    #       })
    #     }, mc.cores = n_cores)
    #   }
    # }
    # 
    # else {
    #   sim_mat <- lapply(1:ncol(classification), function(i) {
    #     sapply(1:ncol(classification), function(j) {
    #       adjustedRandIndex(classification[[i]], classification[[j]])
    #     })
    #   })
    #   
    # }
    # 
    # # Combine the list into a matrix
    # sim_matrix <- do.call(cbind, sim_mat)
    # 
    # # Convert to dissimilarity matrix
    # dissim_matrix <- 1 - sim_matrix
    # dist_mat <- as.dist(dissim_matrix)
    
    #### INSERT IN HERE MDS TEST CODE ####

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
          moc <- parLapply(cl, 1:length(data), function(i) {
            data[[i]] <- as.data.frame(sapply(data[[i]], as.factor))
            m <- do.call(cbind, lapply(data[[i]], function(x) model.matrix(~ x - 1)))
            m <- as.matrix(m)
            colnames(m) <- NULL
            return(m)
          })
          
          # End the cluster
          stopCluster(cl)
        
        } else {
          moc <- mclapply(1:length(data), function(i) {
            data[[i]] <- as.data.frame(sapply(data[[i]], as.factor))
            m <- do.call(cbind, lapply(data[[i]], function(x) model.matrix(~ x - 1)))
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
        for (i in 1:length(data)) {
          data[[i]] <- as.data.frame(sapply(data[[i]], as.factor))
          moc[[i]] <- do.call(cbind, lapply(data[[i]], function(x) model.matrix(~ x - 1)))
          moc[[i]] <- as.matrix(moc[[i]])
          colnames(moc[[i]]) <- NULL
          names(moc)[[i]] <- paste0("MOC - Group", i)  
        }
      }
    }
    return(moc)
  }