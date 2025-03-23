## MultiCOCA
# Load packages 
library(mclust)

# Initialize a list for GMM results
mclust1 <- vector("list", ncol(sim_data))

# Create an empty data frame with the same dimensions as sim_data
sim_data2 <- as.data.frame(matrix(NA, nrow = nrow(sim_data), ncol = ncol(sim_data)))
colnames(sim_data2) <- colnames(sim_data)  # Preserve column names

# Loop through each protein (column) and fit Mclust
for (i in 1:ncol(sim_data)) {  
  mclust1[[i]] <- Mclust(sim_data[, i])
  sim_data2[, i] <- mclust1[[i]]$data  # Store clustered data column-wise
}

