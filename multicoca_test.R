## Full Test - MultiCOCA
# Load packages
library(mclust)
library(NbClust)
library(tidyverse)
library(gtools)

# Load external functions
source("simulateGMM.R")
source("GMMclassifier.R")
source("optimal_hclust.R")
source("constructMOC.R")
source("clusterofclusters.R")
source("multicoca.R")

# Simulate data 
# Initialise parameters
n_groups <- 2
n_clust <- 3
N_col <- 25
n_indiv <- 448
seed <- 4881

# Cluster specific params 
params1 <- list(
  cluster1 = list(mean = rnorm(N_col, mean = -10, sd = 1), cov = cov(matrix(rnorm(N_col*N_col, mean = 0, sd = 0.75), nrow = N_col, ncol = N_col))),
  cluster2 = list(mean = rnorm(N_col, mean = 0, sd = 2), cov = cov(matrix(rnorm(N_col*N_col, mean = 1, sd = 0.35), nrow = N_col, ncol = N_col))),
  cluster3 = list(mean = rnorm(N_col, mean = 10, sd = 0.5), cov = cov(matrix(rnorm(N_col*N_col, mean = -1, sd = 0.15), nrow = N_col, ncol = N_col)))
)

# Generate data
sim_data <- simulateGMM(3, 2, params1, n_indiv = 419, n_col = N_col, 
                          random_seed = 4881, equal_clust = FALSE, equal_groups = FALSE)

data <- sim_data[[1]]
cluster <- sim_data[[2]]

# Classify data using a GMM 
classification <- GMMclassifier(data)

### FROM HERE ONWARDS WILL BE ABSORBED INTO constructMOC FUNCTION

## Calculate dissimilarity matrix
# Construct similarity matrix 
sim_matrix <- matrix(0, nrow = ncol(classification), ncol = ncol(classification),
                     dimnames = list(colnames(classification), colnames(classification)))

set.seed(seed)
for (i in 1:ncol(classification)){
  for(j in 1:ncol(classification)){
    sim_matrix[i,j] <- adjustedRandIndex(classification[[i]], classification[[j]]) # we have called this M (capital)
  }
}

# Pheatmap
pheatmap::pheatmap(sim_matrix)

# Convert to dissimilarity matrix
dissim_matrix <- 1 - sim_matrix
dist_mat <- as.dist(dissim_matrix)

### MDS thresholding would be done here ###

# Figure out best method and linkage for hclust
methods <- list("single","complete","average","ward.D2","ward.D","mcquitty","median","centroid")
indices <- list("cindex","silhouette","dunn","mcclain")

opt_hclust <- optimal_hclust(dist_mat,
                             n_groups,
                             method = methods,
                             index = indices)

# Only frey, mcclain, cindex, sihouette and dunn can be computed. 
# To compute the other indices, data matrix is needed.

# Choose best index and linkage function
opt_hclust$`Number of Clusters`
opt_hclust$`Index Values` # %>% 
  #filter(Index == "dunn" & NClust == ?)

# Apply based on the results above
hclust <- hclust(dist_mat, method = "complete")
opt <- NbClust(diss = dist_mat, distance = NULL, method = "single", 
               min.nc = 2, max.nc = 5, index = "dunn")

# Choose optimal number of clusters
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

### THIS IS WHERE CODE GOING INTO constructMOC ENDS ###

# Create a MOC for each group
moc <- constructMOC(data_group)

# Try MultiCOCA
results <- multicoca(moc, ccClMethod = "hclust", hclustMethod = "complete", 
                        random_seed = 4881, N = 1000, max.iter = 1000, 
                        parallel = TRUE)

# Test against true values
#ARI_group <- groupARI(results, true_clusters)
#print(ARI_group)

# Paul's Code
## THIS FAILS IF MULTICOCA SUSPECTS MORE GROUPS! NEED MORE ROBUST WAY TO CHOOSE NUMBER OF GROUPS
allPermutationsOfTheGroupLabels <- permutations(opt_ng, opt_ng)


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