# Load packages
library(tidyverse)

# Define parameter space
N_col <- c(10, 100, 500, 1000)
param_index <- 1:5
n_groups <- 1:4
models <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI",
            "EEE", "VEE", "EVE", "VVE", "EEV", "VEV", "EVV", "VVV",
            "OTRIMLE")

# Create the grid
grid <- expand.grid(N_col = N_col,
                    param_index = param_index,
                    n_groups = n_groups,
                    models = models,
                    KEEP.OUT.ATTRS = FALSE,
                    stringsAsFactors = FALSE)

# Add a unique task ID (optional, for tracking)
grid$task_id <- seq_len(nrow(grid))

# Subset failed ones 
grid <- grid %>%
  filter(task_id %in% c(25, 105, 143, 185, 269, 305, 325, 445, 485, 505, 529, 
                        549, 586, 605, 618, 645, 665, 805, 829, 845, 849, 897, 
                        909, 925, 949, 969, 978, 998, 1005, 1018, 1038, 1089, 1109))

# Write to CSV
write.table(grid, "input_grid.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 