## Full Test - MultiCOCA
# Load packages
library(lcmm)
library(NbClust)
library(tidyverse)
library(gtools)

# Load external functions
source("simulateLCMM.R")
source("GMMclassifier.R")
source("optimal_hclust.R")
source("constructMOC.R")
source("clusterofclusters.R")
source("multicoca.R")
source("groupARI.R")

# Simulate data 
# Initialise parameters
n_groups <- 3
n_clust <- 3
N_col <- 25
n_indiv <- 448
timepoints <- c(12, 20, 28, 36)
seed <- 4881

# Cluster specific params 
params <- list(
  cluster1 = list(
    fixed_intercept = rnorm(N_col, 8, 1),
    fixed_slope = runif(N_col, -0.1, 0.1),
    random_intercept = runif(N_col, 0.1, 0.5),
    random_slope = runif(N_col, 0.02, 0.15),
    random_cov = matrix(c(mean(runif(N_col, 0.1, 0.5))^2, 0, 0, mean(runif(N_col, 0.02, 0.15))^2), 
                        nrow = 2, byrow = TRUE),
    resid_sd = runif(N_col, 0.1, 0.3)
  ),
  cluster2 = list(
    fixed_intercept = rnorm(N_col, 0, 0.75),
    fixed_slope = runif(N_col, -0.1, 0.1),
    random_intercept = runif(N_col, 0.15, 0.6),
    random_slope = runif(N_col, 0.02, 0.08),
    random_cov = matrix(c(mean(runif(N_col, 0.15, 0.6))^2, 0, 0, mean(runif(N_col, 0.02, 0.08))^2), 
                        nrow = 2, byrow = TRUE),
    resid_sd = runif(N_col, 0.15, 0.25)
  ),
  cluster3 = list(
    fixed_intercept = rnorm(N_col, -8, 0.6),
    fixed_slope = runif(N_col, -0.1, 0.1),
    random_intercept = runif(N_col, 0.05, 0.4),
    random_slope = runif(N_col, 0.04, 0.09),
    random_cov = matrix(c(mean(runif(N_col, 0.05, 0.4))^2, 0, 0, mean(runif(N_col, 0.04, 0.09))^2), 
                        nrow = 2, byrow = TRUE),
    resid_sd = runif(N_col, 0.21, 0.27)
  )
)

# Generate data
sim_data <- simulateLCMM(subject_data = NULL, timepoints, n_clust, n_groups, params, 
                         n_indiv, n_col = N_col, seed, missing = TRUE, timepoint_perc = 0.2, cluster_labels = NULL, 
                         equal_clust = FALSE)

# Seperate out data
group_clusterID <- sim_data[[2]]
group <- sim_data[[3]]
data <- sim_data[[1]]

# Classify data using an LCMM
# classification <- list()
# columns <- names(data %>% dplyr::select(-Subject_ID, -Time))
# max_clusters <- floor(sqrt(length(columns)))
# if (max_clusters < 2) max_clusters <- 2  # at least 2 clusters
# 
# for (col in columns) {
#   numeric_data <- data %>% dplyr::select(Subject_ID, Time, all_of(col))
#   f <- as.formula(paste(col, "~ Time"))
#   
#   # Fit 1-class model
#   model_prev <- lcmm::hlme(
#     fixed = f,
#     random = ~ Time,
#     subject = "Subject_ID",
#     ng = 1,
#     data = numeric_data,
#     nwg = FALSE
#   )
#   
#   bics <- data.frame(Classes = 1, BIC = model_prev$BIC)
#   
#   for (k in 2:max_clusters) {
#       model_k <- lcmm::hlme(
#         fixed = f,
#         random = ~ Time,
#         subject = "Subject_ID",
#         ng = k,
#         mixture = ~ Time,
#         data = numeric_data,
#         nwg = TRUE,
#         B = model_prev
#       )
#     
#     bics <- rbind(bics, data.frame(Classes = k, BIC = model_k$BIC))
#     model_prev <- model_k
#   }
#   
#   best_k <- bics$Classes[which.min(bics$BIC)]
#   message("Best number of classes for ", col, ": ", best_k)
#   
#   classification[[col]] <- list(best_k = best_k, bics = bics)
# }
# classification <- as.data.frame(classification)

col <- columns[1]
numeric_data <- data %>% dplyr::select(Subject_ID, Time, all_of(col))
f <- as.formula(paste(col, "~ Time"))

model1 <- lcmm::hlme(
  fixed = f,
  random = ~ Time,
  subject = "Subject_ID",
  ng = 1,
  data = numeric_data,
  nwg = FALSE
)
print(model1$BIC)

for (i in 2:10) {
model <- lcmm::hlme(
  fixed = f,
  random = ~ Time,
  subject = "Subject_ID",
  ng = i,
  mixture = ~ Time,
  data = numeric_data,
  nwg = TRUE,
  B = model1
)
print(model$BIC)
}

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