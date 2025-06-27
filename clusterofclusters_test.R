## COCA test area
# ==== Setup ====
# Load functions
source("simulateGMM.R")
source("constructMOC.R")
source("clusterofclusters.R")

# Load packages 
library(mclust)

# === Simulation 1: single clustering structure and clearly defined data ====
# Define parameters
N_col <- 10
means <- list(cluster1 = -10, cluster2 = 0, cluster3 = 10)
params1 <- lapply(means, function(m) {
  list(mean = rnorm(N_col, mean = m, sd = 0.5),
       sd = rep(1, N_col))
}) 

# Simulate raw data
sim_data <- simulateGMM(3, 1, params1, n_indiv = 419, n_col = N_col, random_seed = 4881)
clusters1 <- sim_data$`Cluster ID`
sim_data <- sim_data$`Simulated Data`[, 1:N_col]

# Clustering 
classification <- list()
for (i in seq_len(ncol(sim_data))) {
  classification[[i]] <- Mclust(sim_data[,i])$classification
}
classification <- as.data.frame(classification)
names(classification) <- NULL

# Build MOC
moc1 <- constructMOC(classification)

# COCA 
outputCOCA <- clusterofclusters(moc1, N = 10000, max.iter = 10000, parallel_process = TRUE)

# ARI Assessment
clusters1a <- outputCOCA$clusterLabels
adjustedRandIndex(clusters1, clusters1a)

# ==== Simulation 2: single clustering structure and less defined data ====
# Define parameters
N_col <- 10
means <- list(cluster1 = -3, cluster2 = 0, cluster3 = 3)
params2 <- lapply(means, function(m) {
  list(mean = rnorm(N_col, mean = m, sd = 0.5),
       sd = rep(1, N_col))
}) 

# Simulate raw data
sim_data2 <- simulateGMM(3, 1, params2, n_indiv = 419, n_col = N_col, random_seed = 4881)
clusters2 <- sim_data2$`Cluster ID`
sim_data2 <- sim_data2$`Simulated Data`[, 1:N_col]

# Clustering 
classification2 <- list()
for (i in seq_len(ncol(sim_data2))) {
  classification2[[i]] <- Mclust(sim_data2[,i])$classification
}
classification2 <- as.data.frame(classification2)
names(classification2) <- NULL

# Build MOC
moc2 <- constructMOC(classification2)

# COCA 
outputCOCA2 <- clusterofclusters(moc2, N = 1000, max.iter = 1000, parallel_process = TRUE)

# COCA manual
N <- dim(moc2)[1]
maxK <- 10
cm <- array(NA, c(N, N, maxK - 1))
clLabels <- array(NA, c(maxK - 1, N))
for (i in 2:maxK){
  cm[, , i-1] <- coca::consensusCluster(data = moc2, K = i, B = 10000)
  distances <- as.dist(1 - cm[, ,i-1])
  hclust <- hclust(distances, method = "complete")
  clLabels[i-1,] <- cutree(hclust, i)
}
K <- coca::maximiseSilhouette(cm, clLabels, maxK)$K
cm_opt <- cm[, ,K-1]
distances_opt <- as.dist(1 - cm_opt)
print(as.matrix(distances_opt)[1:5, 1:5])
clLabels_opt <- clLabels[K-1,]

# ARI Assessment
clusters2a <- outputCOCA2$clusterLabels
adjustedRandIndex(clusters2, clusters2a)

# ==== Simulation 3: Multiple clustering structures and clearly defined data ====
# Parameters
N_col <- 10
means <- list(cluster1 = -10, cluster2 = 0, cluster3 = 10)
params3 <- lapply(means, function(m) {
  list(mean = rnorm(N_col, mean = m, sd = 0.5),
       sd = rep(1, N_col))
}) 

# Simulate raw data 
sim_data3 <- simulateGMM(3, 2, params1, n_indiv = 419, n_col = N_col, random_seed = 4881)
clusters3 <- sim_data3$`Cluster ID per individual per group`
group3 <- sim_data3$`Group ID`
sim_data3 <- sim_data3$`Simulated Data`[, 1:N_col]

# Clustering 
classification3 <- list()
for (i in seq_len(ncol(sim_data2))) {
  classification3[[i]] <- Mclust(sim_data3[,i])$classification
}
classification3 <- as.data.frame(classification3)
names(classification3) <- NULL

# Build MOC
moc3 <- constructMOC(classification3)

# COCA 
outputCOCA3 <- clusterofclusters(moc3, N = 10000, max.iter = 10000, parallel_process = TRUE)

# ARI Assessment
clusters3a <- outputCOCA3$clusterLabels
adjustedRandIndex(clusters3[[1]], clusters3a)
adjustedRandIndex(clusters3[[2]], clusters3a)
