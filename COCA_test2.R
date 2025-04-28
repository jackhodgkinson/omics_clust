## COCA test 
# Load packages
library(mclust)
library(NbClust)

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
# clusters <- sample(c(1,2,3), 450, replace = TRUE, prob = {p <- runif(length(c(1,2,3))); p / sum(p)})
# data1 <- simulateGMM(3, 2, params1, n_indiv = 419, n_col = N_col,
#                      random_seed = seed, cluster_labels = clusters,
#                      equal_clust = FALSE, equal_groups = FALSE)
# data1 <- data1[[1]]
data1 <- simulateGMM(3, 2, params1, n_indiv = 419, n_col = N_col,
                     random_seed = seed,
                     equal_clust = FALSE, equal_groups = FALSE)
true_clusters <- data1[[2]]
true_groups <- data1[[3]]
data1 <- data1[[1]]
# N_col <- 10
# params1 <- list(
#   cluster1 = list(mean = runif(N_col, -10, -5), sd = runif(N_col, 0.5, 1)),
#   cluster2 = list(mean = runif(N_col, 10, 15), sd = runif(N_col, 0.75, 3)),
#   cluster3 = list(mean = rep(0, N_col), sd = rep(5, N_col)))
# data2 <- simulateGMM(3, 1, params1, n_indiv = 450, n_col = N_col,
#                      random_seed = 9377, cluster_labels = clusters,
#                      equal_clust = FALSE, equal_groups = FALSE)
# data2 <- data2[[1]]
# data <- list(data1, data2)
data <- data1

# Data from COCA package
# seed <- 4881
# data <- list()
# data[[1]] <- as.matrix(read.csv(system.file("extdata", "dataset1.csv",
#                                             package = "coca"), row.names = 1))
# data[[2]] <- as.matrix(read.csv(system.file("extdata", "dataset2.csv",
#                                             package = "coca"), row.names = 1))
# data[[3]] <- as.matrix(read.csv(system.file("extdata", "dataset3.csv",
#                                             package = "coca"), row.names = 1))
# 
# true_labels <- as.matrix(read.csv(system.file("extdata", "cluster_labels.csv",
#                                               package = "coca"), row.names = 1))
# 
# N_col <- ncol(data[[1]])
# N_row <- nrow(data[[1]])

# Cluster each column if single dataset
if (class(data) != "list") {
  
  # Create empty list for mclust and classification data frame 
  mclust <- vector("list", ncol(data))
  classification <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))
                                  
  # Convert data frame to a list of columns 
  data2 <- as.list(data)
  
  # Fit Mclust to each protein to obtain classification
  set.seed(seed)
  
  for (i in 1:length(data2)) {
    colnames(data2[[i]]) <- NULL
    mclust[[i]] <- Mclust(data2[[i]])
    classification[, i] <- mclust[[i]]$classification  # Store clustered data column-wise
    }

# Cluster each column if a list of datasets
} else {
    
  # Initialize lists to store Mclust results and classification matrices for each data frame
  mclust_results <- list()  
  classification_results <- list() 
  
  for (i in 1:length(data)) {
    
    # Assign data to dataframe i 
    data2 <- data[[i]]
    
    # Create empty list for mclust and classification data frame 
    mclust <- vector("list", ncol(data2))
    classification <- as.data.frame(matrix(NA, nrow = nrow(data2), ncol = ncol(data2)))
                                    
    # Convert data frame to a list of columns 
    data3 <- as.list(data2)
                                    
    # Fit Mclust to each protein to obtain classification 
    set.seed(seed)
    for (j in 1:length(data3)) {
      colnames(data3[[j]]) <- NULL
      mclust[[j]] <- Mclust(data3[[j]])
      classification[, j] <- mclust[[j]]$classification  # Store clustered data column-wise
    }
    mclust_results[[paste0("Dataset",i)]] <- mclust
    classification_results[[paste0("Dataset",i)]] <- classification
  }
}
                              
# When do we combine the different sources of information? When do we need to 
# only keep the rows that share ID var? 

# We need to think about how to handle missing data - could just remove if < 10% 
# and stop the function if larger than this and advise user to investigate.


# Construct similarity matrix 
sim_matrix <- matrix(0, nrow = ncol(classification), ncol = ncol(classification),
                     dimnames = list(colnames(classification), colnames(classification)))

set.seed(seed)
for (i in 1:ncol(classification)){
  for(j in 1:ncol(classification)){
    sim_matrix[i,j] <- adjustedRandIndex(classification[[i]], classification[[j]])
  }
}

## Pheatmap 
pheatmap::pheatmap(sim_matrix)

# Convert to dissimilarity matrix
dissim_matrix <- 1 - sim_matrix
dist_mat <- as.dist(dissim_matrix)

# Here we should add a threshold whether the data needs splitting into groups or not! 
# Need a criterion to decide whether there is more than one cluster or not. 

# Apply hclust to dissim_matrix - complete gives same as intended input. 
hclust <- hclust(dist_mat) # Try different linkage options - single, complete, average, ward.D2

# Plot hclust
plot(hclust)

# Assess impact in simulation study. 

# Figure out optimal number of groups 
set.seed(seed)
opt <- NbClust(diss = dist_mat, distance = NULL, method = "single", 
               min.nc = 2, max.nc = 9, index = "dunn")

# Try different linkage functions with different indices and see which is most accurate. Report this in thesis. 
if (is.vector(opt$Best.nc)) {
  opt_ng <- as.numeric(opt$Best.nc[1])
} else {
  opt_ng <- as.numeric(names(which.max(table(opt$Best.nc[1, ]))))
}

# Split dendogram
groups <- cutree(hclust, k = opt_ng)

# Test grouping 
if (exists("true_groups")) {
  adjustedRandIndex(groups, true_groups)
}

# Split the columns into datasets by group 
data_group <- list()

for(i in 1:length(unique(groups))){
  data_group[[i]] <- classification[groups == i]
  names(data_group)[[i]] <- paste0("Group",i)
}

# Create a MOC for each data group
moc <- constructMOC(data_group)

# Pheatmap 
pheatmap::pheatmap(moc$`MOC - Group1`)
pheatmap::pheatmap(moc$`MOC - Group2`)

# For COCA function, how do we ensure that the max k is reasonable? Same as hclust problem.

# Feed into COCA for each data group 
results <- multicoca(moc, random_seed = 4881, N = 1000, max.iter = 1000)

# Test against true values
ari_g1 <- adjustedRandIndex(true_clusters$group1_clusterid, results$Group1$clusterLabels)
ari_g2 <- adjustedRandIndex(true_clusters$group2_clusterid, results$Group2$clusterLabels)

# Figure out why this is not working!