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

## Test: single clustering structure
example1 <- simulateGMM(3, 1, params1, n_indiv = 419, n_col = N_col, random_seed = 4881)

# Seperate data into data and clusters
clusters <- example1$cluster
example1 <- example1[, 1:N_col]

# Assign column names
colnames(example1) <- paste0("Protein", 1:N_col)

# Assign labels to the clusters
clusters <- factor(clusters, labels = c("1","2","Unassigned"))

# Produce heatmap 
example1_new <- as.matrix(example1)
annotationRow <- as.data.frame(clusters)
names(annotationRow) <- "Clusters"
rownames(example1_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
pheatmap::pheatmap(example1_new, annotation_row = annotationRow)

## Test: two clustering structures defined randomly 
example2 <- simulateGMM(3, 2, params1, n_indiv = 419, n_col = N_col, random_seed = 4881)

# Seperate data into data and clusters
clusters <- example2$cluster
example2 <- example2[, 1:N_col]

# Assign column names
colnames(example2) <- paste0("Protein", 1:N_col)

# Assign labels to the clusters
clusters <- factor(clusters, labels = c("1","2","Unassigned"))

# Produce heatmap 
example2_new <- as.matrix(example2)
annotationRow <- as.data.frame(clusters)
names(annotationRow) <- "Clusters"
rownames(example2_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
pheatmap::pheatmap(example2_new, annotation_row = annotationRow)

## Test: two clustering structures defined randomly 
example3 <- simulateGMM(3, 2, params1, n_indiv = 419, n_col = N_col, random_seed = 4881, grouping = "hclust")

# Seperate data into data and clusters
clusters <- example3$cluster
example3 <- example3[, 1:N_col]

# Assign column names
colnames(example3) <- paste0("Protein", 1:N_col)

# Assign labels to the clusters
clusters <- factor(clusters, labels = c("1","2","Unassigned"))

# Produce heatmap 
example3_new <- as.matrix(example3)
annotationRow <- as.data.frame(clusters)
names(annotationRow) <- "Clusters"
rownames(example3_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
pheatmap::pheatmap(example3_new, annotation_row = annotationRow)

## Parameters - means close together
N_col <- 10
params2 <- list(
  cluster1 = list(mean = runif(N_col, -20, 10), sd = runif(N_col, 0.5, 5)),  
  cluster2 = list(mean = runif(N_col, -20, 10), sd = runif(N_col, 1, 4)), 
  cluster3 = list(mean = rep(0, N_col), sd = rep(5, N_col)))  

## Test: single clustering structure
example1a <- simulateGMM(3, 1, params2, n_indiv = 419, n_col = N_col, random_seed = 4881)

# Seperate data into data and clusters
clusters1a <- example1a$cluster
example1a <- example1a[, 1:N_col]

# Assign column names
colnames(example1a) <- paste0("Protein", 1:N_col)

# Assign labels to the clusters
clusters1a <- factor(clusters1a, labels = c("1","2","Unassigned"))

# Produce heatmap 
example1a_new <- as.matrix(example1a)
annotationRow <- as.data.frame(clusters1a)
names(annotationRow) <- "Clusters"
rownames(example1a_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
pheatmap::pheatmap(example1a_new, annotation_row = annotationRow)

## Test: two clustering structures defined randomly 
example2a <- simulateGMM(3, 2, params2, n_indiv = 419, n_col = N_col, random_seed = 4881)

# Seperate data into data and clusters
clusters2a <- example2a$cluster
example2a <- example2a[, 1:N_col]

# Plot the data and colour by cluster for the first 5 clusters
example2a_data_plot <- example2a %>%
  select(Protein1:Protein5) %>%
  cbind(clusters2a)

ggpairs(example2a_data_plot, aes(colour = clusters2a))

# Produce heatmap 
example2a_new <- as.matrix(example2a)
annotationRow <- as.data.frame(clusters2a)
names(annotationRow) <- "Clusters"
rownames(example2a_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
pheatmap::pheatmap(example2a_new, annotation_row = annotationRow)

## Test: two clustering structures defined randomly 
example3a <- simulateGMM(3, 2, params2, n_indiv = 419, n_col = N_col, random_seed = 4881, grouping = "hclust")

# Seperate data into data and clusters
clusters3a <- example3a$cluster
example3a <- example3a[, 1:N_col]

# Assign column names
colnames(example3a) <- paste0("Protein", 1:N_col)

# Assign labels to the clusters
clusters3a <- factor(clusters3a, labels = c("1","2","Unassigned"))

# Produce heatmap 
example3a_new <- as.matrix(example3a)
annotationRow <- as.data.frame(clusters3a)
names(annotationRow) <- "Clusters"
rownames(example3a_new) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
pheatmap::pheatmap(example3a_new, annotation_row = annotationRow)

