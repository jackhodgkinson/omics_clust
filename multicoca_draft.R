## MultiCOCA
# Load packages 
library(mclust)
library(coca)

# MultiCOCA function 
multicoca <- function(data, clusters) {

# Initialize a list for GMM results
mclust1 <- vector("list", ncol(data))

# Create an empty data frame with the same dimensions as sim_data
data2 <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))

# Loop through each protein (column) and fit Mclust
for (i in 1:ncol(data)) {  
  mclust1[[i]] <- Mclust(data[, i])
  data2[, i] <- mclust1[[i]]$classification  # Store clustered data column-wise
}

# How do we deal if Mclust comes up with different numbers of clusters? 

# Change data to MOC
#moc <- coca::buildMOC(data2, M = ncol(data2), K = clusters)

# Try to create MOC using function and then try to make MOC myself to avoid double clustering.
# Visualise the MOC
# 

return(data2)
}

data2 <- multicoca(sim_data, 2)



