#### SIMULATION STUDY ####

### SETUP
## Load packages
library(coca)
library(MASS)
library(parallel)
library(gtools)

## Load functions
source("GMMclassifier.R")
source("numCores.R")
source("simulateGMM.R")

# Command line argument for array jobs
args <- commandArgs(trailingOnly = TRUE)
N_col <- as.numeric(args[1])
param_index <- as.numeric(args[2])
n_clust <- as.numeric(args[3])
n_groups <- as.numeric(args[4])
linkage <- args[5]
index <- args[6]
task_id <- as.numeric(args[7])
seed <- task_id
set.seed(task_id)

# Create many different parameters 
params_list <- list(
  # 1. Well-separated means, identity covariance
  list(
    cluster1 = list(mean = rnorm(N_col, mean = -3, sd = 0.1), cov = diag(N_col)),
    cluster2 = list(mean = rnorm(N_col, mean = 0,  sd = 0.1), cov = diag(N_col)),
    cluster3 = list(mean = rnorm(N_col, mean = 3,  sd = 0.1), cov = diag(N_col))
  ),
  
  # 2. Closer means, identity covariance
  list(
    cluster1 = list(mean = rnorm(N_col, mean = -1.5, sd = 0.1), cov = diag(N_col)),
    cluster2 = list(mean = rnorm(N_col, mean = 0,    sd = 0.1), cov = diag(N_col)),
    cluster3 = list(mean = rnorm(N_col, mean = 1.5,  sd = 0.1), cov = diag(N_col))
  ),
  
  # 3. Closer still, diagonal covariance with varied scale
  list(
    cluster1 = list(mean = rnorm(N_col, mean = -1, sd = 0.1), cov = diag(runif(N_col, 0.5, 1.5))),
    cluster2 = list(mean = rnorm(N_col, mean = 0,  sd = 0.1), cov = diag(runif(N_col, 0.5, 1.5))),
    cluster3 = list(mean = rnorm(N_col, mean = 1,  sd = 0.1), cov = diag(runif(N_col, 0.5, 1.5)))
  ),
  
  # 4. Full covariance, small overlap
  list(
    cluster1 = list(mean = rnorm(N_col, mean = -0.75, sd = 0.1), cov = cov(matrix(rnorm(N_col*N_col, 0, 0.5), nrow = N_col, ncol = N_col))),
    cluster2 = list(mean = rnorm(N_col, mean = 0,     sd = 0.1), cov = cov(matrix(rnorm(N_col*N_col, 0, 0.5), nrow = N_col, ncol = N_col))),
    cluster3 = list(mean = rnorm(N_col, mean = 0.75,  sd = 0.1), cov = cov(matrix(rnorm(N_col*N_col, 0, 0.5), nrow = N_col, ncol = N_col)))
  ),
  
  # 5. High overlap, full covariance + correlation
  list(
    cluster1 = list(mean = rnorm(N_col, mean = -0.3, sd = 0.1), cov = cov(matrix(rnorm(N_col*N_col, 1, 0.3), nrow = N_col, ncol = N_col))),
    cluster2 = list(mean = rnorm(N_col, mean = 0,    sd = 0.1), cov = cov(matrix(rnorm(N_col*N_col, 1, 0.3), nrow = N_col, ncol = N_col))),
    cluster3 = list(mean = rnorm(N_col, mean = 0.3,  sd = 0.1), cov = cov(matrix(rnorm(N_col*N_col, 1, 0.3), nrow = N_col, ncol = N_col)))
  )
)

param_labels <- c(
  "Well-separated means, identity cov",
  "Closer means, identity cov",
  "Closer still, diag cov varied scale",
  "Full covariance, small overlap",
  "High overlap, full covariance + correlation"
)

# Generate parameters 
param_set <- params_list[[param_index]] 
param_label <- param_labels[[param_index]]

# Get cores
n_cores <- numCores()

# Generate empty results data frame
results <- data.frame(
  Parameter.ID = NULL,
  Parameter.Label = NULL,
  Columns = NULL,
  `True Number of Views` = NULL,
  Linkage = NULL,
  Index = NULL,
  `Determined Number of Views` = NULL,
  `ARI of Views` = NULL, 
  `ARI of Individuals per View` = NULL,
  RunTime = NULL,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

### STAGE 1: GENERATE DATA 
# Generate data from a GMM
data <- simulateGMM(n_clust, n_groups, param_set, n_indiv = 419, n_col = N_col,
                    random_seed = seed,
                    equal_clust = FALSE, equal_groups = FALSE)

# Obtain the true cluster labels
true_clusters <- data[[2]]

# Obtain the true view labels
if (n_groups > 1) {
  true_groups <- data[[3]]
}

# Obtain the true data
data <- data[[1]]

### STAGE 2: CLUSTERING OF INDIVIDUALS
# Fit MClust to the individuals per column / dataset.
classification <- GMMclassifier(data)

### STAGE 3: CLUSTERING OF VIEWS
# 


### STAGE 4: CLUSTER OF CLUSTERS ANALYSIS 


### STAGE 5: ASSESSMENT OF CLUSTERING
## Permutation of Group Labels
opt_ng <- allPermutationsOfTheGroupLabels <- permutations(opt_ng, opt_ng)


ariResultsMatrix <- matrix(nrow = nrow(allPermutationsOfTheGroupLabels), ncol = opt_ng)
for(i in 1:nrow(allPermutationsOfTheGroupLabels))
{
  currentPermutation       <- allPermutationsOfTheGroupLabels[i,]
  currentOrderingOfResults <- results[currentPermutation]
  
  for(j in 1:opt_ng){
    ariResultsMatrix[i,j] <- adjustedRandIndex(currentOrderingOfResults[[j]]$clusterLabels, true_clusters[,j])
  }
  
}

print(ariResultsMatrix)

bestMatch <- which.max(rowSums(ariResultsMatrix))

finalARI  <- ariResultsMatrix[bestMatch,]

print("ARI for the derived clustering by group vs true clustering by group")
print(finalARI)





