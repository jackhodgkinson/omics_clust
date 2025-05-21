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

# Write to CSV
write.table(grid, "input_grid.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)