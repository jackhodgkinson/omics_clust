# Load packages 
library(tidyverse)
library(GGally)

# Set parameters
set.seed(4881)
n <- 419
p <- 10

# Simulate partiicpant cluster membership 
indiv_clust <- sample(c(1, 2, 3), size = n, replace = TRUE)
table(indiv_clust)

# Create matrix to store data
sim_data <- matrix(NA, ncol = p, nrow = n)

# Generate random distribution parameters per cluster
cluster_params <- list(
  cluster_1 = list(mean = runif(p, -20, 10), sd = runif(p, 0.5, 5)),  
  cluster_2 = list(mean = runif(p, -20, 10), sd = runif(p, 1, 4)), 
  cluster_3 = list(mean = rep(0, p), sd = rep(5, p))   
)

# Generate data
for (protein in 1:p) {
  
  protein_data <- numeric(n)
  
  for (i in 1:n) {
    
    # Simulate random normals for each protein and cluster
    mu <- cluster_params[[paste0("cluster_", indiv_clust[i])]]$mean[[protein]]
    sigma <-  cluster_params[[paste0("cluster_", indiv_clust[i])]]$sd[[protein]]
    
    # Generate cluster assignment 
    protein_data[i] <- rnorm(1, mean = mu, sd = sigma)
  }
  
  # Add protein data to new column
  sim_data[, protein] <- protein_data
}

# Convert to data frame 
sim_data <- as.data.frame(sim_data)

# Assign column names
colnames(sim_data) <- paste0("Protein", 1:p)

# Assign labels to the clusters
cluster <- factor(indiv_clust, labels = c("1","2","Unassigned"))

# Plot the data and colour by cluster for the first 5 clusters
sim_data_plot <- sim_data %>%
  select(Protein1:Protein5) 

ggpairs(sim_data_plot, aes(color = cluster))
