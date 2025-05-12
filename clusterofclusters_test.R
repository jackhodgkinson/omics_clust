## COCA test area
# Load functions
source("simulateGMM.R")
source("constructMOC.R")
source("clusterofclusters.R")

# Load packages 
library(mclust)

# Create data - single clustering structure and clearly defined data
N_col <- 10
means <- list(cluster1 = -10, cluster2 = 0, cluster3 = 10)
params1 <- lapply(means, function(m) {
  list(mean = rnorm(N_col, mean = m, sd = 0.5),
       sd = rep(1, N_col))
}) 
sim_data <- simulateGMM(3, 1, params1, n_indiv = 419, n_col = N_col, random_seed = 4881)
clusters1 <- sim_data$`Cluster ID`
sim_data <- sim_data$`Simulated Data`[, 1:N_col]

outputCOCA <- clusterofclusters(sim_data, max.iter = 10000)
clusters2 <- outputCOCA$clusterLabels
adjustedRandIndex(clusters1, clusters2)

# Create data - single clustering structure and less defined data
N_col <- 10
params2 <- list(
  cluster1 = list(mean = rnorm(N_col, mean = -2, sd = 1), sd = rep(2, N_col)),
  cluster2 = list(mean = rnorm(N_col, mean = 0, sd = 0.5), sd = rep(2.2, N_col)),
  cluster3 = list(mean = rnorm(N_col, mean = 2, sd = 0.6), sd = rep(1.5, N_col))
)
sim_data2 <- simulateGMM(3, 1, params2, n_indiv = 419, n_col = N_col, random_seed = 4881)
clusters1a <- sim_data2$`Cluster ID`
sim_data2 <- sim_data2$`Simulated Data`[, 1:N_col]

outputCOCA_2 <- clusterofclusters(sim_data2, max.iter = 10000)
clusters2a <- outputCOCA$clusterLabels
adjustedRandIndex(clusters1a, clusters2a)

# Create data - multiple clustering structures and clearly defined data

# sim_data2 <- simulateGMM(3, 2, params1, n_indiv = 419, n_col = N_col, random_seed = 4881)
# clusters1a <- sim_data2$`Cluster ID per individual per group`
# group1a <- sim_data2$`Group ID`
# sim_data2 <- sim_data2$`Simulated Data`[, 1:N_col]
# 
# moc <- constructMOC(sim_data2)
# outputCOCA_2 <- clusterofclusters(sim_data2)
# clusters2a <- outputCOCA_2$clusterLabels
# adjustedRandIndex(clusters1a, clusters2a)
