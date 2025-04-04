# Test examples from the COCA package (Cabassi and Kirk (2019))
# Turn on warning messages
options(warn = 1)

# Load function
source("clusterofclusters.R")
source("simulateGMM.R")

# Load package 
library(coca)

# Load data
data <- list()
data[[1]] <- as.matrix(read.csv(system.file("extdata", "dataset1.csv",
                                            package = "coca"), row.names = 1))
data[[2]] <- as.matrix(read.csv(system.file("extdata", "dataset2.csv",
                                            package = "coca"), row.names = 1))
data[[3]] <- as.matrix(read.csv(system.file("extdata", "dataset3.csv",
                                            package = "coca"), row.names = 1))
true_labels <- as.matrix(read.csv(system.file("extdata", "cluster_labels.csv",
                                              package = "coca"), row.names = 1))

# Run original COCA function 
outputBuildMOC <- coca::buildMOC(data, M = 3, K = 5, distances = "cor")
moc <- outputBuildMOC$moc
example <- coca::coca(moc, K = 5)
pred_labels <- example$clusterLabels

# Run my COCA function 
example_new <- clusterofclusters(data, k = 5)
pred_labels_new <- example$clusterLabels

# Compare results 
(ari = adjustedRandIndex(true_labels, pred_labels))
(ari_new = adjustedRandIndex(true_labels, pred_labels_new))

# Try original coca on my simulated data 
N_col <- 10
params1 <- list(
  cluster1 = list(mean = runif(N_col, -20, -10), sd = runif(N_col, 0.5, 1)),  
  cluster2 = list(mean = runif(N_col, 10, 20), sd = runif(N_col, 0.5, 2)), 
  cluster3 = list(mean = runif(N_col, 0, 0.5), sd = runif(N_col, 1, 10)))  
sim_data <- simulateGMM(3, 1, params1, n_indiv = 419, n_col = N_col, random_seed = 4881, grouping = "random")
clusters1 <- sim_data$cluster
sim_data <- sim_data[, 1:N_col]
sim_list <- lapply(sim_data, function(x) data.frame(x))

outputBuildMOC_sim <- buildMOC(sim_list, M = N_col, distances = "euclidean")
moc_sim <- outputBuildMOC_sim$moc
outputCOCA_sim <- coca::coca(moc)
clusterLabels <- outputCOCA_sim$clusterLabels