# Simulate longitudinal data
# ==== Setup ====
# Load packages
library(mclust)
library(tidyverse)
library(parallel)

# Load external functions
source("simulateLCMM.R")
source("numCores.R")

# Simulated data 
seed <- 4881
set.seed(seed)
N_col <- c(10,100,500,1000)
n_groups <- c(2, 3, 4, 5)
n_clust <- c(2, 3, 4, 5)
n_indiv <- 325

# Parallel 
parallel_process <- TRUE

# Specify parameters 
params <- list(
  cluster1 = list(mean = rnorm(N_col, mean = -1.5, sd = 0.1), cov = diag(N_col)),
  cluster2 = list(mean = rnorm(N_col, mean = 0,    sd = 0.1), cov = diag(N_col)),
  cluster3 = list(mean = rnorm(N_col, mean = 1.5,  sd = 0.1), cov = diag(N_col))
)

# ==== Simulation ====
# Set seed
set.seed(seed)

# Simulation loop
sim_total <- 50
datasets_per_sim <- 0  # Will be set after first sim

# Log file
log_file <- "data_sim_long.txt"
cat("Simulation of Longitudinal Data \n", file = log_file)

# Create data simulation folder
output_dir <- "simulated_data"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Dataset and seed counter
data_counter <- 0
seed_counter <- 0

# Simulation loop
for (sim in 1:sim_total) {
  
  for (n in seq(N_col)) {
    
    if (N_col[n] == 10) {
      n_groups <- c(2, 3)
      n_clust <- c(2, 3)
    } else if (N_col[n] == 100) {
      n_groups <- c(2, 3, 4)
      n_clust <- c(2, 3, 4)
    } else {
      n_groups <- c(2, 3, 4, 5)
      n_clust <- c(2, 3, 4, 5)
    }
    
    for (g in seq(n_groups)) {
      for (c in seq(n_clust)) {
        
        # Generate unique seed based on base + counter
        data_counter <- data_counter + 1
        this_seed <- seed + seed_counter
        set.seed(this_seed)
        
        log_entry <- paste0("data", counter, ".csv \nSimulation: ", sim, ", Seed: ", this_seed, "\n", "\n")
        cat(log_entry, file = log_file, append = TRUE)
        
        # Your simulation using `this_seed`
        
        
        # Save data
        filename <- paste0("data", counter, ".csv")
        full_path <- file.path(output_dir, filename)
        write_csv(sim_data, full_path, row.names = FALSE)
        
        # Reset counter
        seed_counter <- seed_counter + 1
      }
    }
  }
}