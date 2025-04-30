## COCA test area
# Load functions
source("simulateGMM.R")
source("constructMOC.R")
source("clusterofclusters.R")

# Load packages 
library(mclust)

# Create data - single clustering structure 
N_col <- 10
params1 <- list(
  cluster1 = list(mean = runif(N_col, -20, -10), sd = runif(N_col, 0.5, 1)),  
  cluster2 = list(mean = runif(N_col, 10, 20), sd = runif(N_col, 1, 5)), 
  cluster3 = list(mean = runif(N_col, 0, 0.5), sd = runif(N_col, 1, 10)))  
sim_data <- simulateGMM(3, 1, params1, n_indiv = 419, n_col = N_col, random_seed = 4881)
clusters1 <- sim_data$`Cluster ID`
sim_data <- sim_data$`Simulated Data`[, 1:N_col]

# Test on data1
outputCOCA <- clusterofclusters(sim_data, k = 3, max.iter = 10000)
clusters2 <- outputCOCA$clusterLabels
adjustedRandIndex(clusters1, clusters2)

# Create data - multiple clustering structures
sim_data2 <- simulateGMM(3, 2, params1, n_indiv = 419, n_col = N_col, random_seed = 4881)
clusters1a <- sim_data2$`Cluster ID per individual per group`
group1a <- sim_data2$`Group ID`
sim_data2 <- sim_data2$`Simulated Data`[, 1:N_col]

# Test on simulated data 
moc <- constructMOC(sim_data2)
outputCOCA_2 <- clusterofclusters(sim_data2)
clusters2a <- outputCOCA_2$clusterLabels
adjustedRandIndex(clusters1a, clusters2a)
