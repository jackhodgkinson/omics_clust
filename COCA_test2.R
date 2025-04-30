## COCA test 
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
# clusters <- sample(c(1,2,3), 450, replace = TRUE, prob = {p <- runif(length(c(1,2,3))); p / sum(p)})
# data1 <- simulateGMM(3, 2, params1, n_indiv = 419, n_col = N_col,
#                      random_seed = seed, cluster_labels = clusters,
#                      equal_clust = FALSE, equal_groups = FALSE)
# data1 <- data1[[1]]
n_groups <- 2
data1 <- simulateGMM(3, n_groups, params1, n_indiv = 419, n_col = N_col,
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

# Apply hclust to dissim_matrix - try out different linkage functions 
method <- list("single","complete","average","ward.D2","ward.D","mcquitty","median","centroid")

# Set up plot layout 
par(mfrow = c(3, 3), mar = c(3, 3, 3, 2))

# Loop over each method 
for (m in method){
  try({
    hc <- hclust(dist_mat, method = m)
    plot(hc, main = paste("Method:", m), xlab = "", sub = "")
  }, silent = TRUE)
}

# Reset plotting layout
par(mfrow = c(1, 1))

# Do best hclust 
hclust <- hclust(dist_mat)

# Figure out optimal number of groups 
set.seed(seed)
index <- list("cindex","silhouette","dunn","mcclain")

# Try different linkage functions with different indices and see which is most accurate. Report this in thesis.
results <- data.frame(Method = character(), 
                      Index = character(), 
                      NGroups = integer(), 
                      stringsAsFactors = FALSE)

all_index_list <- list()

for (m in method) {
  for (i in index) {
    try({
      opt <- NbClust(diss = dist_mat, distance = NULL,
                     method = m, min.nc = 2, max.nc = 9,
                     index = i)
      best_nc <- opt$Best.nc[[1]]
      opt_results <- rbind(results, data.frame(Method = m,
                                           Index = i,
                                           BestNC = best_nc))
      index_df <- data.frame(NClust = as.numeric(names(opt$All.index)),
                             Value = as.numeric(opt$All.index),
                             Method = m,
                             Index = i)
      all_index_list[[paste(m, i, sep = "_")]] <- index_df
    }, silent = TRUE)
  }
}

# Join with true number of groups 
opt_results <- cbind(opt_results, n_groups)

# Create data frame for plotting
all_index_df <- bind_rows(all_index_list)

# Plot graph 
ggplot(all_index_df, aes(x = NClust, y = Value)) + 
  geom_line() + 
  facet_grid(Method ~ Index, scales = "free_y") + 
  labs(title = "NbClust Index Scores for Each Method and Index",
       x = "Number of Groups",
       y = "Index Value")

# Subset for dunn 
dunn_vals <- all_index_df %>%
  filter(Index == "dunn" & NClust == 2)

# Based off results, calculate optimal number of groups of clusters
opt <- NbClust(diss = dist_mat, distance = NULL, method = "single", 
               min.nc = 2, max.nc = 9, index = "dunn")


# Choose optimal 
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

# Do I need to recluster in groups here again? Or use the existing clusters?
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
results1 <- clusterofclusters(moc$`MOC - Group1`, ccClMethod = "hclust", hclustMethod = "complete", random_seed = 4881, N = 1000, max.iter = 1000)
results2 <- clusterofclusters(moc$`MOC - Group1`, ccClMethod = "kmeans", random_seed = 4881, N = 10000, max.iter = 10000)

# Test COCA against true 
ari_g1a <- adjustedRandIndex(true_clusters$group1_clusterid, results1$clusterLabels)
ari_g1b <- adjustedRandIndex(true_clusters$group1_clusterid, results2$clusterLabels) # Method and number of iterations makes no difference! 

pheatmap::pheatmap(results1$consensusMatrix)
pheatmap::pheatmap(results2$consensusMatrix)

results_hc <- multicoca(moc, ccClMethod = "hclust", hclustMethod = "complete", 
                        random_seed = 4881, N = 1000, max.iter = 1000, 
                        parallel = FALSE) #This gives group 2 k = 2
results_hc_p <- multicoca(moc, ccClMethod = "hclust", hclustMethod = "complete", 
                        random_seed = 4881, N = 1000, max.iter = 1000, ccDistHC = "pearson",
                        parallel = TRUE) #This gives group 2 k = 2
results_hc_s <- multicoca(moc, ccClMethod = "hclust", hclustMethod = "complete", 
                        random_seed = 4881, N = 1000, max.iter = 1000, ccDistHC = "spearman",
                        parallel = TRUE) #This gives group 2 k = 2
results_km <- multicoca(moc, ccClMethod = "kmeans", random_seed = 4881, 
                        N = 2000, max.iter = 2000) # Gives group 2 k = 3

pheatmap::pheatmap(results_hc$Group1$consensusMatrix)
pheatmap::pheatmap(results_hc$Group2$consensusMatrix)
pheatmap::pheatmap(results_km$Group1$consensusMatrix)
pheatmap::pheatmap(results_km$Group2$consensusMatrix)

# Test against true values
ari_g1_hc <- adjustedRandIndex(true_clusters$group1_clusterid, results_hc$Group1$clusterLabels)
ari_g2_hc <- adjustedRandIndex(true_clusters$group2_clusterid, results_hc$Group2$clusterLabels)

ari_g1_hc_p <- adjustedRandIndex(true_clusters$group1_clusterid, results_hc_p$Group1$clusterLabels)
ari_g2_hc_p <- adjustedRandIndex(true_clusters$group2_clusterid, results_hc_p$Group2$clusterLabels)

ari_g1_hc_s <- adjustedRandIndex(true_clusters$group1_clusterid, results_hc_s$Group1$clusterLabels)
ari_g2_hc_s <- adjustedRandIndex(true_clusters$group2_clusterid, results_hc_s$Group2$clusterLabels)

ari_g1_km <- adjustedRandIndex(true_clusters$group1_clusterid, results_km$Group1$clusterLabels)
ari_g2_km <- adjustedRandIndex(true_clusters$group2_clusterid, results_km$Group2$clusterLabels)

# Spearman correlation and Euclidean distance are better than random guessing for group 2! 

# Figure out why this is not working!