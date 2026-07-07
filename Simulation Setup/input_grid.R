# Load packages
library(tidyverse)

# Define other parameter values
N_col <- c(10, 100, 500, 1000)
param_index <- 1:5
methods <- c("single", "complete", "average", "ward.D2", "ward.D", "mcquitty", "median", "centroid")
indices <- c("cindex", "silhouette", "dunn", "mcclain")

# Function to determine n_clust based on N_col
get_cluster_range <- function(ncol) {
  if (ncol <= 10) {
    return(2:3)
  } else if (ncol <= 100) {
    return(2:4)
  } else {
    return(2:5)
  }
}

# Function to determine n_groups based on N_col
get_group_range <- function(ncol) {
  if (ncol <= 10) {
    return(2)
  } else if (ncol <= 100) {
    return(2:3)
  } else {
    return(2:5)
  }
}

# Build grid dynamically
grid <- map_dfr(N_col, function(ncol) {
  n_clust <- get_cluster_range(ncol)
  n_groups <- get_group_range(ncol)
  
  expand.grid(
    N_col = ncol,
    param_index = param_index,
    n_clust = n_clust,
    n_groups = n_groups,
    methods = methods,
    indices = indices,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
})

# Add task ID
grid <- grid %>%
  mutate(task_id = row_number())

# Write to CSV
write.table(grid, "input_grid.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)