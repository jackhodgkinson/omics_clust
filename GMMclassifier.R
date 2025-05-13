# GMMclassifier.R
# Load numCores.R function
source("numCores.R")

# Function
GMMclassifier <- function(data,
                          parallel = TRUE) {
  
  # Detect number of cores
  n_cores <- numCores()
  
  # Cluster each column if single dataset
  if (class(data) != "list") {
    
    # Create empty list for mclust and classification data frame 
    mclust <- vector("list", ncol(data))
    classification <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))
    
    # Convert data frame to a list of columns 
    data2 <- as.list(data)
    
    # Set seed
    set.seed(seed)
    
    if (parallel) {
    
      if (.Platform$OS.type == "windows"){
      
        # Set up parallel cluster
        cl <- makeCluster(n_cores)
        clusterExport(cl, ls(envir = environment()), envir = environment())
        
        # Fit Mclust to each protein to obtain classification
        classification_results <- parLapply(cl, data2, function(i) {
          colnames(i) <- NULL
          return(Mclust(i)$classification)
        })
        
        # Combine into a classification data frame
        classification <- as.data.frame(do.call(cbind, classification_results))
      
        # End parallel cluster
        stopCluster(cl) 
        
    } else {
        
        # Fit Mclust to each protein to obtain classification
        classification_results <- mclapply(data2, function(i) {
          colnames(i) <- NULL
          return(mclust::Mclust(i)$classification)
        },
        mc.cores = n_cores)
        
        # Combine into a classification data frame
        classification <- as.data.frame(do.call(cbind, classification_results))
        
    }
    } else {
          
          # Fit Mclust to each protein to obtain classification
          classification_results <- lapply(data2, function(i) {
            colnames(i) <- NULL
            return(Mclust(i)$classification)
          })
          
          # Combine into a classification data frame
          classification <- as.data.frame(do.call(cbind, classification_results))
        }
      }
    
    # Cluster each column if a list of datasets
    else {
    
    # Initialize lists to store Mclust results and classification matrices for each data frame
    mclust_results <- list()  
    classification_results <- list() 
    
    # Set seed 
    set.seed(seed)
    
    if (parallel) {
    
      if (.Platform$OS.type == "windows"){
      
        # Set parallel cluster
        clo <- makeCluster(n_cores)
        clusterExport(clo, ls(envir = environment()), envir = environment())
        
        # Loop over each element in data list
        data_process <- function(data) {
          
          # Create empty list for mclust and classification data frame 
          mclust <- vector("list", ncol(data))
          classification <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))
          
          # Convert data frame to a list of columns 
          data2 <- as.list(data)
          
          # Set up inner parallel cluster for columns
          set.seed(seed)
          cli <- makeCluster(n_cores)
          clusterExport(cli, ls(envir = environment()), envir = environment())
          clusterEvalQ(cli, library(mclust))
          
          # Fit Mclust to each protein to obtain classification 
          class_res_in <- parLapply(cli, data2, function(i) {
            colnames(i) <- NULL
            mclust <- Mclust(i)
            return(mclust$classification)  # Store clustered data column-wise
          })
          
          # Stop the cluster
          stopCluster(cli)
          
          # Combine classification results into a data frame
          classification <- as.data.frame(do.call(cbind, class_res_in))
          return(classification)
        }
        
        # Parallelize the processing of each dataset in the 'data' list (outer parallelization)
        classification_results <- parLapply(clo, data, function(dataset) {
          process_dataset(dataset)  # Process each dataset in parallel
        })
        
        # Stop the outer cluster after the work is done
        stopCluster(clo)
        
        # Name the results for each dataset
        names(classification_results) <- paste0("Dataset", 1:length(data))
        
    } else {
      
      # Loop over each element in data list
      data_process <- function(data) {
        
        # Create empty list for mclust and classification data frame 
        mclust <- vector("list", ncol(data))
        classification <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))
        
        # Convert data frame to a list of columns 
        data2 <- as.list(data)
        
        # Set seed
        set.seed(seed)
        
        # Fit Mclust to each protein to obtain classification 
        class_res_in <- mclapply(data2, function(i) {
          colnames(i) <- NULL
          mclust <- Mclust(i)
          return(mclust$classification)},
          mc.cores = n_cores)
        
        # Combine classification results into a data frame
        classification <- as.data.frame(do.call(cbind, class_res_in))
        return(classification)
      }
      
      # Parallelize the processing of each dataset in the 'data' list (outer parallelization)
      classification_results <- mclapply(data, function(dataset) {process_dataset(dataset)},
                                         mc.cores = n_cores)
      
      # Name the results for each dataset
      names(classification_results) <- paste0("Dataset", 1:length(data))
        
    }
      
    } else {
    
      # Loop over each element in data list
      data_process <- function(data) {
        
        # Create empty list for mclust and classification data frame 
        mclust <- vector("list", ncol(data))
        classification <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))
        
        # Convert data frame to a list of columns 
        data2 <- as.list(data)
        
        # Fit Mclust to each protein to obtain classification 
        class_res_in <- lapply(data2, function(i) {
          colnames(i) <- NULL
          mclust <- Mclust(i)
          return(mclust$classification)  # Store clustered data column-wise
        })
        
        # Combine classification results into a data frame
        classification <- as.data.frame(do.call(cbind, class_res_in))
        return(classification)
      }
      
      # Parallelize the processing of each dataset in the 'data' list (outer parallelization)
      classification_results <-lapply(data, function(dataset) {
        process_dataset(dataset)  # Process each dataset in parallel
      })
      
      # Name the results for each dataset
      names(classification_results) <- paste0("Dataset", 1:length(data))
      classification <- classification_results
      
    }
    }
  
  return(classification)
}