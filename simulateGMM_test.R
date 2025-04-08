#### Testing of simulateGMM function in different scenarios 
### Install relevant packages 
library(tidyverse)
library(GGally)
library(mclust)

### Load function 
source("simulateGMM.R")

### Testing
## Parameters - means far apart 
N_col <- 10
params1 <- list(
  cluster1 = list(mean = runif(N_col, -10, -5), sd = runif(N_col, 0.5, 1)),  
  cluster2 = list(mean = runif(N_col, 10, 15), sd = runif(N_col, 0.75, 3)), 
  cluster3 = list(mean = rep(0, N_col), sd = rep(5, N_col)))  

## Test: single clustering structure with equal sized clusters
example1 <- simulateGMM(3, 1, params1, n_indiv = 419, n_col = N_col, random_seed = 4881)

# Seperate data into data and clusters
clusters <- example1[[2]]
example1_data <- example1[[1]]

# Assign column names
colnames(example1_data) <- paste0("Protein", 1:N_col)

# Assign labels to the clusters
clusters <- factor(clusters, labels = c("1","2","Unassigned"))

# Produce heatmap 
example1_new <- as.matrix(example1_data)
annotationRow <- as.data.frame(clusters)
names(annotationRow) <- "Clusters"
rownames(example1_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
pheatmap::pheatmap(example1_new, annotation_row = annotationRow)

## Test: single clustering structure with different sized clusters 
example1.1 <- simulateGMM(3, 1, params1, n_indiv = 419, n_col = N_col, 
                        random_seed = 4881, equal_clust = FALSE)

# Seperate data into data and clusters
clusters <- example1.1[[2]]
example1.1_data <- example1.1[[1]]

# Assign column names
colnames(example1.1_data) <- paste0("Protein", 1:N_col)

# Assign labels to the clusters
clusters <- factor(clusters, labels = c("1","2","Unassigned"))

# Produce heatmap 
example1.1_new <- as.matrix(example1.1_data)
annotationRow <- as.data.frame(clusters)
names(annotationRow) <- "Clusters"
rownames(example1.1_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
pheatmap::pheatmap(example1.1_new, annotation_row = annotationRow)

## Test: two clustering structures with equal sized clusters and groups 
example2 <- simulateGMM(3, 2, params1, n_indiv = 419, n_col = N_col, random_seed = 4881)

# Seperate data into data and clusters
example2_data <- example2[[1]]
clusters <- example2[[2]]

# Assign column names
colnames(example2_data) <- paste0("Protein", 1:N_col)

# Assign labels to the clusters
clusters[] <- lapply(clusters, function(x) factor(x, labels = c("1","2","Unassigned")))

# Produce heatmap - focused on group 1
example2_new <- as.matrix(example2_data)
annotationRow <- as.data.frame(clusters[1])
names(annotationRow) <- "Clusters"

annotationCol <- as.data.frame(as.factor(example2[[3]]))
names(annotationCol) <- "ProteinClusters"
rownames(annotationCol) <- colnames(example2_new)

rownames(example2_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
annotationCol$ProteinClusters <- as.factor(annotationCol$ProteinCluster)
pheatmap::pheatmap(example2_new[order(clusters[[1]]),
                                order(annotationCol$ProteinClusters)], 
                   annotation_row = annotationRow,
                   annotation_col = annotationCol[order(annotationCol$ProteinClusters), , drop = FALSE], 
                   cluster_rows = F)

# Produce heatmap - focused on group 2
example2_new <- as.matrix(example2_data)
annotationRow <- as.data.frame(clusters[2])
names(annotationRow) <- "Clusters"

annotationCol <- as.data.frame(as.factor(example2[[3]]))
names(annotationCol) <- "ProteinClusters"
rownames(annotationCol) <- colnames(example2_new)

rownames(example2_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
annotationCol$ProteinClusters <- as.factor(annotationCol$ProteinClusters)
pheatmap::pheatmap(example2_new[order(clusters[[2]]),
                                order(annotationCol$ProteinClusters)], 
                   annotation_row = annotationRow,
                   annotation_col = annotationCol[order(annotationCol$ProteinClusters), , drop = FALSE], 
                   cluster_rows = F)

## Test: two clustering structures with different sized clusters and groups 
example2.1 <- simulateGMM(3, 2, params1, n_indiv = 419, n_col = N_col, 
                        random_seed = 4881, equal_clust = FALSE, equal_groups = FALSE)

# Seperate data into data and clusters
example2.1_data <- example2.1[[1]]
clusters <- example2.1[[2]]

# Assign column names
colnames(example2.1_data) <- paste0("Protein", 1:N_col)

# Assign labels to the clusters
clusters[] <- lapply(clusters, function(x) factor(x, labels = c("1","2","Unassigned")))

# Produce heatmap - focused on group 1
example2.1_new <- as.matrix(example2.1_data)
annotationRow <- as.data.frame(clusters[1])
names(annotationRow) <- "Clusters"

annotationCol <- as.data.frame(example2.1[[3]])
names(annotationCol) <- "ProteinClusters"
rownames(annotationCol) <- colnames(example2.1_new)

rownames(example2.1_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
annotationCol$ProteinClusters <- as.factor(annotationCol$ProteinClusters)
pheatmap::pheatmap(example2.1_new[order(clusters[["group1_clusterid"]]),
                                  order(annotationCol$ProteinClusters)], 
                   annotation_row = annotationRow,
                   annotation_col = annotationCol[order(annotationCol$ProteinClusters), , drop = FALSE], 
                   cluster_rows = F,
                   cluster_cols = F)

# Produce heatmap - focused on group 2
example2.1_new <- as.matrix(example2.1_data)
annotationRow <- as.data.frame(clusters[2])
names(annotationRow) <- "Clusters"

annotationCol <- as.data.frame(example2.1[[3]])
names(annotationCol) <- "ProteinClusters"
annotationCol$ProteinClusters <- factor(annotationCol$ProteinClusters)
rownames(annotationCol) <- colnames(example2.1_new)

rownames(example2.1_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
pheatmap::pheatmap(example2.1_new[order(clusters[["group2_clusterid"]]),
                                  order(annotationCol$ProteinClusters)], 
                   annotation_row = annotationRow,
                   annotation_col = annotationCol[order(annotationCol$ProteinClusters), , drop = FALSE], 
                   cluster_rows = F,
                   cluster_cols = F)


## Test: three clustering structures with equal sized clusters and groups  
example3 <- simulateGMM(3, 3, params1, n_indiv = 419, n_col = N_col, random_seed = 4881)

# Seperate data into data and clusters
example3_data <- example3[[1]]
clusters <- example3[[2]]

# Assign column names
colnames(example3_data) <- paste0("Protein", 1:N_col)

# Assign labels to the clusters
clusters[] <- lapply(clusters, function(x) factor(x, labels = c("1","2","Unassigned")))

# Produce heatmap - focused on group 1
example3_new <- as.matrix(example3_data)
annotationRow <- as.data.frame(clusters[1])
names(annotationRow) <- "Clusters"

annotationCol <- as.data.frame(as.factor(example3[[3]]))
names(annotationCol) <- "ProteinClusters"
rownames(annotationCol) <- colnames(example3_new)

rownames(example3_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
annotationCol$ProteinClusters <- as.factor(annotationCol$ProteinClusters)
pheatmap::pheatmap(example3_new[order(clusters[["group1_clusterid"]]),
                                order(annotationCol$ProteinClusters)], 
                   annotation_row = annotationRow,
                   annotation_col = annotationCol[order(annotationCol$ProteinClusters), , drop = FALSE], 
                   cluster_rows = F)

# Produce heatmap - focused on group 2
example3_new <- as.matrix(example3_data)
annotationRow <- as.data.frame(clusters$group2_clusterid)
names(annotationRow) <- "Clusters"

annotationCol <- as.data.frame(as.factor(example3[[3]]))
names(annotationCol) <- "ProteinClusters"
rownames(annotationCol) <- colnames(example3_new)

rownames(example3_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
annotationCol$ProteinClusters <- as.factor(annotationCol$ProteinClusters)
pheatmap::pheatmap(example3_new[order(clusters[["group2_clusterid"]]),
                                order(annotationCol$ProteinClusters)], 
                   annotation_row = annotationRow,
                   annotation_col = annotationCol[order(annotationCol$ProteinClusters), , drop = FALSE], 
                   cluster_rows = F)

# Produce heatmap - focused on group 3
example3_new <- as.matrix(example3_data)
annotationRow <- as.data.frame(clusters$group3_clusterid)
names(annotationRow) <- "Clusters"

annotationCol <- as.data.frame(as.factor(example3[[3]]))
names(annotationCol) <- "ProteinClusters"
rownames(annotationCol) <- colnames(example3_new)

rownames(example3_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
annotationCol$ProteinClusters <- as.factor(annotationCol$ProteinClusters)
pheatmap::pheatmap(example3_new[order(clusters[["group3_clusterid"]]),
                                order(annotationCol$ProteinClusters)], 
                   annotation_row = annotationRow,
                   annotation_col = annotationCol[order(annotationCol$ProteinClusters), , drop = FALSE], 
                   cluster_rows = F)

## Test: three clustering structures with different sized clusters and groups  
example3.1 <- simulateGMM(3, 3, params1, n_indiv = 419, n_col = N_col, 
                          random_seed = 4881, equal_clust = FALSE, equal_groups = FALSE)

# Seperate data into data and clusters
example3.1_data <- example3.1[[1]]
clusters <- example3.1[[2]]

# Assign column names
colnames(example3.1_data) <- paste0("Protein", 1:N_col)

# Assign labels to the clusters
clusters[] <- lapply(clusters, function(x) factor(x, labels = c("1","2","Unassigned")))

# Produce heatmap - focused on group 1
example3.1_new <- as.matrix(example3.1_data)
annotationRow <- as.data.frame(clusters[1])
names(annotationRow) <- "Clusters"

annotationCol <- as.data.frame(as.factor(example3.1[[3]]))
names(annotationCol) <- "ProteinClusters"
rownames(annotationCol) <- colnames(example3.1_new)

rownames(example3.1_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
annotationCol$ProteinClusters <- as.factor(annotationCol$ProteinClusters)
pheatmap::pheatmap(example3.1_new[order(clusters[["group1_clusterid"]]),
                                order(annotationCol$ProteinClusters)], 
                   annotation_row = annotationRow,
                   annotation_col = annotationCol[order(annotationCol$ProteinClusters), , drop = FALSE], 
                   cluster_rows = F,
                   cluster_cols = F)

# Produce heatmap - focused on group 2
example3.1_new <- as.matrix(example3.1_data)
annotationRow <- as.data.frame(clusters$group2_clusterid)
names(annotationRow) <- "Clusters"

annotationCol <- as.data.frame(as.factor(example3.1[[3]]))
names(annotationCol) <- "ProteinClusters"
rownames(annotationCol) <- colnames(example3.1_new)

rownames(example3.1_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
annotationCol$ProteinClusters <- as.factor(annotationCol$ProteinClusters)
pheatmap::pheatmap(example3.1_new[order(clusters[["group2_clusterid"]]),
                                order(annotationCol$ProteinClusters)], 
                   annotation_row = annotationRow,
                   annotation_col = annotationCol[order(annotationCol$ProteinClusters), , drop = FALSE], 
                   cluster_rows = F,
                   cluster_cols = F)

# Produce heatmap - focused on group 3
example3.1_new <- as.matrix(example3.1_data)
annotationRow <- as.data.frame(clusters$group3_clusterid)
names(annotationRow) <- "Clusters"

annotationCol <- as.data.frame(as.factor(example3.1[[3]]))
names(annotationCol) <- "ProteinClusters"
rownames(annotationCol) <- colnames(example3_new)

rownames(example3.1_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
annotationCol$ProteinClusters <- as.factor(annotationCol$ProteinClusters)
pheatmap::pheatmap(example3.1_new[order(clusters[["group3_clusterid"]]),
                                order(annotationCol$ProteinClusters)], 
                   annotation_row = annotationRow,
                   annotation_col = annotationCol[order(annotationCol$ProteinClusters), , drop = FALSE], 
                   cluster_rows = F,
                   cluster_cols = F)