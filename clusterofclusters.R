## clusterofclusters.R
clusterofclusters <- function(data,              # Input as data frame
                              N = 2000,          # Number of iterations of Consensus Clustering step
                              max.iter = 2000    # Maximum number of iterations for k-means clustering
                              ) 
  {
  
  # Load packages - write into dependencies when writing package.
  library(mclust)
  library(coca)
  
  # Initialize a list for model-based clustering results and number of clusters
  mclust1 <- vector("list", ncol(data))
  G <- vector("numeric", ncol(data))
  
  # Create an empty data frame with the same dimensions as sim_data
  classification <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))
  
  # Loop through each column and fit Mclust
  for (i in 1:ncol(data)) {
    mclust1[[i]] <- Mclust(data[, i])
    classification[, i] <- mclust1[[i]]$classification  # Store clustered data column-wise
  }
  
  # Convert to factor (to ensure proper binary encoding)
  classification <- as.data.frame(lapply(classification, as.factor))
  
  # Create the MOC
  moc <- do.call(cbind, lapply(classification, function(x) model.matrix(~ x - 1)))
  moc <- as.matrix(moc)
  colnames(moc) <- NULL
  
  # Fit Cluster of Cluster Analysis
  outputCOCA <- coca::coca(moc, B = N, maxIterKM = max.iter)

  return(outputCOCA)
}