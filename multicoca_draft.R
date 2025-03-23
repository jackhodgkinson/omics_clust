## MultiCOCA
# Load packages 
library(mclust)

# MultiCOCA function 
multicoca <- function(data) {

# Initialize a list for GMM results
mclust1 <- vector("list", ncol(data))

# Create an empty data frame with the same dimensions as sim_data
data2 <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))
colnames(data2) <- colnames(data)  # Preserve column names

# Loop through each protein (column) and fit Mclust
for (i in 1:ncol(data)) {  
  mclust1[[i]] <- Mclust(data[, i])
  data2[, i] <- mclust1[[i]]$data  # Store clustered data column-wise
}

}


