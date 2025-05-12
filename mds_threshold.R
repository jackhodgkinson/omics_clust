## MDS Thresholding 
# Turn off graphics
graphics.off()

# Load packages
library(mclust)
library(NbClust)
library(tidyverse)
library(parallel)

# Load external functions
source("simulateGMM.R")
source("constructMOC.R")
source("clusterofclusters.R")
source("multicoca.R")

# Simulated data 
seed <- 4881
set.seed(seed)
N_col <- 10

# Define parameters for data simulation
# Well defined clusters
params1 <- list(
  cluster1 = list(mean = rnorm(N_col, mean = -10, sd = 0.6), sd = rep(0.25, N_col)),
  cluster2 = list(mean = rnorm(N_col, mean = 0, sd = 0.4), sd = rep(0.3, N_col)),
  cluster3 = list(mean = rnorm(N_col, mean = 10, sd = 0.7), sd = rep(0.15, N_col))
  )

# Less defined clusters with correlation
params2 <- list(
  cluster1 = list(mean = rnorm(N_col, mean = -4, sd = 0.6), cov = cov(matrix(rnorm(N_col*N_col, mean = 0, sd = 0.5), nrow = N_col, ncol = N_col))),
  cluster2 = list(mean = rnorm(N_col, mean = 0, sd = 0.4), cov = cov(matrix(rnorm(N_col*N_col, mean = 1, sd = 0.25), nrow = N_col, ncol = N_col))),
  cluster3 = list(mean = rnorm(N_col, mean = 4, sd = 0.7), cov = cov(matrix(rnorm(N_col*N_col, mean = -1, sd = 0.15), nrow = N_col, ncol = N_col)))
)

# N group of clusters
n_groups <- 3
data <- simulateGMM(3, n_groups, params2, n_indiv = 419, n_col = N_col,
                     random_seed = seed,
                     equal_clust = FALSE, equal_groups = FALSE)
true_clusters <- data[[2]]
true_groups <- data[[3]]
data <- data[[1]]

# Plot data
data <- as.matrix(data)
annotationRow <- as.data.frame(true_clusters[1])
names(annotationRow) <- "Clusters"

annotationCol <- as.data.frame(as.factor(true_groups))
names(annotationCol) <- "ProteinClusters"
rownames(annotationCol) <- colnames(data)

rownames(data) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
annotationCol$ProteinClusters <- as.factor(annotationCol$ProteinCluster)
true_clust <- as.numeric(true_clusters[["group1_clusterid"]])
ordered_data <- data[order(true_clust), order(annotationCol$ProteinClusters)]
pheatmap::pheatmap(
  ordered_data,  
  annotation_row = annotationRow,  
  annotation_col = annotationCol[order(annotationCol$ProteinClusters), , drop = FALSE],  
  cluster_rows = FALSE  
)

# Cluster each column if single dataset
if (class(data) != "list") {
  
  # Create empty list for mclust and classification data frame 
  mclust <- vector("list", ncol(data))
  classification <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))
  
  # Convert data frame to a list of columns 
  data2 <- as.list(data)
  
  # Set seed
  set.seed(seed)
  
  if (.Platform$OS.type == "windows"){
  
    # Set up parallel cluster
    cl <- makeCluster(detectCores() - 1)
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
    mc.cores = detectCores() - 1)
    
    # Combine into a classification data frame
    classification <- as.data.frame(do.call(cbind, classification_results))
  }
  
  # Cluster each column if a list of datasets
} else {
  
  # Initialize lists to store Mclust results and classification matrices for each data frame
  mclust_results <- list()  
  classification_results <- list() 
  
  # Set seed 
  set.seed(seed)
  
  if (.Platform$OS.type == "windows"){
  
    # Set parallel cluster
    clo <- makeCluster(detectCores() - 1)
    clusterExport(clo, ls(envir = environment()), envir = environment())
    clusterEvalQ(clo, library(mclust))
    
    # Loop over each element in data list
    data_process <- function(data) {
      
      # Create empty list for mclust and classification data frame 
      mclust <- vector("list", ncol(data))
      classification <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))
      
      # Convert data frame to a list of columns 
      data2 <- as.list(data)
      
      # Set up inner parallel cluster for columns
      set.seed(seed)
      cli <- makeCluster(detectCores() - 1)
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
      mc.cores = detectCores() - 1)
    
    # Combine classification results into a data frame
    classification <- as.data.frame(do.call(cbind, class_res_in))
    return(classification)
  }
  
  # Parallelize the processing of each dataset in the 'data' list (outer parallelization)
  classification_results <- mclapply(data, function(dataset) {process_dataset(dataset)},
                                     mc.cores = detectCores() - 1)
  
  # Name the results for each dataset
  names(classification_results) <- paste0("Dataset", 1:length(data))
    
  }
}

# Construct empty similarity matrix 
sim_matrix <- matrix(0, nrow = ncol(classification), ncol = ncol(classification),
                     dimnames = list(colnames(classification), colnames(classification)))

# Set seed 
set.seed(seed)

if (.Platform$OS.type == "windows") {
  cl2 <- makeCluster(detectCores() - 1)
  clusterExport(cl2, ls(envir = environment()), envir = environment())
  
  # Parallelised computation of similarity matrix
  sim_mat <- parLapply(cl, 1:ncol(classification), function(i) {
    sapply(1:ncol(classification), function(j) {
      adjustedRandIndex(classification[[i]], classification[[j]])
    })
  })
  
  # Stop the parallel cluster after the work is done
  stopCluster(cl2)
  
} else {
  sim_mat <- mclapply(1:ncol(classification), function(i) {
    sapply(1:ncol(classification), function(j) {
      adjustedRandIndex(classification[[i]], classification[[j]])
    })
  })
}

# Combine the list into a matrix
sim_matrix <- do.call(cbind, sim_mat)

# Pheatmap 
pheatmap::pheatmap(sim_matrix)

# Convert to dissimilarity matrix
dissim_matrix <- 1 - sim_matrix
dist_mat <- as.dist(dissim_matrix)

# Perform MDS on the distance matrix
mds <- cmdscale(dist_mat)

# Fit Gaussian and Gaussian Mixture (2 components) to MDS
model <- Mclust(mds, G=2)$modelName
if (!is.null(model)) {
  test <- mclust::mclustBootstrapLRT(mds, modelName = model, nboot = 1000, level = 0.95, maxG = 1)
  p <- test$p.value
}

# Produce MDS plot 
mds_data <- as.data.frame(mds)

# Seperate if p >= 0.05 
if (p < 0.05){
mds_plot <- ggplot(mds_data, aes(x = V1, y = V2)) + 
  geom_point() + 
  labs(x = "MDS1", y = "MDS2",
       title = "Multidimensional Scaling (MDS) Plot of Distance Matrix")
mds_plot
} else {
  
  
  
  # Once grouped append the group to the mds_data and color by the group
  mds_data$Group <- 
  mds_plot <- ggplot(mds_data, aes(x = V1, y = V2, color = )) + 
    geom_point() + 
    labs(x = "MDS1", y = "MDS2",
         title = "Multidimensional Scaling (MDS) Plot of Distance Matrix")
  mds_plot
}
