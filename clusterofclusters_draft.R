## clusterofclusters.R
clusterofclusters <- function(data,              # Input as data frame or list of data frames
                              k = 1:9,           # Number of clusters
                              N = 1000,          # Number of iterations of Consensus Clustering step
                              max.iter = 1000    # Maximum number of iterations for k-means clustering
                              ) 
  {
  
  # Load packages - write into dependencies when writing package.
  library(mclust)
  library(coca)
  
  # Convert data frame into a list of data frame columns if not already a list of data frames
  if(class(data) != "list"){
    n_col <- ncol(data)
    n_row <- nrow(data)
    data_l <- lapply(data, function(x) data.frame(x))
  }
  else{
    n_row <- max(sapply(data, nrow))
    n_col <- length(data)
    data_l <- data
  }

  # Initialize a list for model-based clustering results and number of clusters
  mclust1 <- vector("list", n_col)
  G <- vector("numeric", n_col)

  # Create an empty data frame with the same dimensions as sim_data
  classification <- as.data.frame(matrix(NA, nrow = n_row, ncol = n_col))

  # Loop through each list element and fit Mclust
  for (i in 1:length(data_l)) {
    colnames(data_l[[i]]) <- NULL
    mclust[[i]] <- Mclust(data_l[[i]], G = k)
    classification[, i] <- mclust1[[i]]$classification  # Store clustered data column-wise
  }

  # Convert to factor (to ensure proper binary encoding)
  classification <- as.data.frame(lapply(classification, as.factor))
  
  # Here I need to work out the clustering of the proteins. (LOOK AT PROTOCOL - HCLUST?)

  # Create the MOC - this needs to be per protein group. single MOC function with an option!
  moc <- do.call(cbind, lapply(classification, function(x) model.matrix(~ x - 1)))
  moc <- as.matrix(moc)
  colnames(moc) <- NULL

  # Fit Cluster of Cluster Analysis
  outputCOCA <- coca::coca(moc, B = N, maxIterKM = max.iter, verbose = FALSE)

  return(outputCOCA)
}