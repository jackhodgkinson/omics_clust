## Testing cophenetic correlation coefficient 
# Delete all plots
graphics.off()

# Load packages
library(mclust)
library(NbClust)
library(tidyverse)

# Load external functions
source("simulateGMM.R")
source("constructMOC.R")
source("clusterofclusters.R")
source("multicoca.R")

# Simulated data 
seed <- 4881
set.seed(seed)
N_col <- 10
params1 <- list(
  cluster1 = list(mean = runif(N_col, -10, -5), sd = runif(N_col, 0.5, 1)),
  cluster2 = list(mean = runif(N_col, 10, 15), sd = runif(N_col, 0.75, 3)),
  cluster3 = list(mean = rep(0, N_col), sd = rep(5, N_col)))


# Test for different number of clusters 
n_groups <- c(1:5)
for (i in n_groups) {
  data <- simulateGMM(3, i, params1, n_indiv = 419, n_col = N_col,
                      random_seed = seed,
                      equal_clust = FALSE, equal_groups = FALSE)
  true_clusters <- data[[2]]
  data <- data[[1]]
  
  # Cluster each column if single dataset
  if (class(data) != "list") {
    
    # Create empty list for mclust and classification data frame 
    mclust <- vector("list", ncol(data))
    classification <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))
    
    # Convert data frame to a list of columns 
    data2 <- as.list(data)
    
    # Fit Mclust to each protein to obtain classification
    set.seed(seed)
    
    for (j in 1:length(data2)) {
      colnames(data2[[j]]) <- NULL
      mclust[[j]] <- Mclust(data2[[j]])
      classification[, j] <- mclust[[j]]$classification  # Store clustered data column-wise
    }
    
    # Cluster each column if a list of datasets
  } else {
    
    # Initialize lists to store Mclust results and classification matrices for each data frame
    mclust_results <- list()  
    classification_results <- list() 
    
    for (k in 1:length(data)) {
      
      # Assign data to dataframe i 
      data2 <- data[[i]]
      
      # Create empty list for mclust and classification data frame 
      mclust <- vector("list", ncol(data2))
      classification <- as.data.frame(matrix(NA, nrow = nrow(data2), ncol = ncol(data2)))
      
      # Convert data frame to a list of columns 
      data3 <- as.list(data2)
      
      # Fit Mclust to each protein to obtain classification 
      set.seed(seed)
      for (l in 1:length(data3)) {
        colnames(data3[[l]]) <- NULL
        mclust[[l]] <- Mclust(data3[[l]])
        classification[, l] <- mclust[[l]]$classification  # Store clustered data column-wise
      }
      mclust_results[[paste0("Dataset",k)]] <- mclust
      classification_results[[paste0("Dataset",k)]] <- classification
    }
  }
  
  # Construct similarity matrix 
  sim_matrix <- matrix(0, nrow = ncol(classification), ncol = ncol(classification),
                       dimnames = list(colnames(classification), colnames(classification)))
  
  set.seed(seed)
  for (m in 1:ncol(classification)){
    for(n in 1:ncol(classification)){
      sim_matrix[m,n] <- adjustedRandIndex(classification[[m]], classification[[n]])
    }
  }
  
  ## Pheatmap 
  pheatmap::pheatmap(sim_matrix)
  
  # Convert to dissimilarity matrix
  dissim_matrix <- 1 - sim_matrix
  dist_mat <- as.dist(dissim_matrix)
  
  # # Apply hclust to dissim_matrix - try out different linkage functions 
  # method <- list("single","complete","average","ward.D2","ward.D","mcquitty","median","centroid")
  # 
  # # Set up plot layout 
  # par(mfrow = c(3, 3), mar = c(3, 3, 3, 2))
  # 
  # # Loop over each method 
  # for (m in method){
  #   try({
  #     hc <- hclust(dist_mat, method = m)
  #     plot(hc, main = paste("Method:", m), xlab = "", sub = "")
  #   }, silent = TRUE)
  # }
  # 
  # # Reset plotting layout
  # par(mfrow = c(1, 1))
  
  # Do MDS 
  fit <- cmdscale(dist_mat, k = 2)
  plot(fit)
  fit_mc <- Mclust(fit, G = 1:10)
  plot(fit_mc$BIC)
  class <- fit_mc$classification
  
}