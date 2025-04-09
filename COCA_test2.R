## COCA test 
# Load packages
library(mclust)
library(NbClust)

# Load external functions
source("simulateGMM.R")
source("constructMOC.R")

# Simulate data 
N_col <- 10
params1 <- list(
  cluster1 = list(mean = runif(N_col, -10, -5), sd = runif(N_col, 0.5, 1)),  
  cluster2 = list(mean = runif(N_col, 10, 15), sd = runif(N_col, 0.75, 3)), 
  cluster3 = list(mean = rep(0, N_col), sd = rep(5, N_col)))  
sim_data <- simulateGMM(3, 2, params1, n_indiv = 419, n_col = N_col, 
                        random_seed = 4881, equal_clust = FALSE, equal_groups = FALSE)
true_clusters <- sim_data[[2]]
true_groups <- sim_data[[3]]
sim_data <- sim_data[[1]]

# Create empty list for mclust and classification data frame 
mclust <- vector("list", N_col)
classification <- as.data.frame(matrix(NA, nrow = 419, ncol = N_col))

# Fit Mclust to each protein to obtain classification 
for (i in 1:length(sim_data)) {
  colnames(sim_data[[i]]) <- NULL
  mclust[[i]] <- Mclust(sim_data[[i]])
  classification[, i] <- mclust[[i]]$classification  # Store clustered data column-wise
}

# Construct similarity matrix 
sim_matrix <- matrix(0, nrow = ncol(classification), ncol = ncol(classification),
                     dimnames = list(colnames(classification), colnames(classification)))

for (i in 1:ncol(classification)){
  for(j in 1:ncol(classification)){
    sim_matrix[i,j] <- adjustedRandIndex(classification[[i]], classification[[j]])
  }
}

# Convert to dissimilarity matrix
dissim_matrix <- 1 - sim_matrix
dist_mat <- as.dist(dissim_matrix)

# Here we should add a threshold whether the data needs splitting into groups or not! 

# Apply hclust to dissim_matrix - complete gives same as intended input. 
hclust <- hclust(dist_mat)

# Figure out optimal number of clusters - "all" seem to be most accurate - "alllong" takes too long. 
opt <- NbClust(classification, distance = "euclidean", method = "complete", 
               min.nc = 2, max.nc = 10, index = "all")
opt_nc <- as.numeric(names(which.max(table(opt$Best.nc[1, ]))))

# Split dendogram
groups <- cutree(hclust, k = opt_nc)

# Test grouping 
adjustedRandIndex(groups, true_groups)

# Split the columns into datasets by group 
data_group <- list()

for(i in 1:length(unique(groups))){
  data_group[[i]] <- classification[groups == i]
  names(data_group)[[i]] <- paste0("Group",i)
}

# Create a MOC for each data group
moc <- constructMOC(data_group)

# Feed into COCA for each data group 
results <- list()
for (i in 1:length(moc)) {
  results[[i]] <- coca::coca(moc[[i]], B = 1000, maxIterKM = 1000, verbose = FALSE)
  names(results)[[i]] <- paste0("Group",i)
}

adjustedRandIndex(true_clusters$group1_clusterid, results$Group1$clusterLabels)
adjustedRandIndex(true_clusters$group2_clusterid, results$Group2$clusterLabels)

# Very high ARI - how do we ensure the ARIs are the same per group? Is that required?  
