## MDS Thresholding 
# Turn off graphics
graphics.off()

# Load packages
library(mclust)
library(NbClust)
library(tidyverse)
library(parallel)
library(otrimle)

# Load external functions
source("simulateGMM.R")
source("GMMclassifier.R")
source("numCores.R")

# Command line argument for array jobs
args <- commandArgs(trailingOnly = TRUE)
N_col <- as.numeric(args[1])
param_index <- as.numeric(args[2])
n_groups <- as.numeric(args[3])
model <- args[4]
print(model)
task_id <- as.numeric(args[5])
seed <- task_id
set.seed(task_id)

# Create many different parameters 
params_list <- list(
  # 1. Well-separated means, identity covariance
  list(
    cluster1 = list(mean = rnorm(N_col, mean = -3, sd = 0.1), cov = diag(N_col)),
    cluster2 = list(mean = rnorm(N_col, mean = 0,  sd = 0.1), cov = diag(N_col)),
    cluster3 = list(mean = rnorm(N_col, mean = 3,  sd = 0.1), cov = diag(N_col))
  ),
  
  # 2. Closer means, identity covariance
  list(
    cluster1 = list(mean = rnorm(N_col, mean = -1.5, sd = 0.1), cov = diag(N_col)),
    cluster2 = list(mean = rnorm(N_col, mean = 0,    sd = 0.1), cov = diag(N_col)),
    cluster3 = list(mean = rnorm(N_col, mean = 1.5,  sd = 0.1), cov = diag(N_col))
  ),
  
  # 3. Closer still, diagonal covariance with varied scale
  list(
    cluster1 = list(mean = rnorm(N_col, mean = -1, sd = 0.1), cov = diag(runif(N_col, 0.5, 1.5))),
    cluster2 = list(mean = rnorm(N_col, mean = 0,  sd = 0.1), cov = diag(runif(N_col, 0.5, 1.5))),
    cluster3 = list(mean = rnorm(N_col, mean = 1,  sd = 0.1), cov = diag(runif(N_col, 0.5, 1.5)))
  ),
  
  # 4. Full covariance, small overlap
  list(
    cluster1 = list(mean = rnorm(N_col, mean = -0.75, sd = 0.1), cov = cov(matrix(rnorm(N_col*N_col, 0, 0.5), nrow = N_col, ncol = N_col))),
    cluster2 = list(mean = rnorm(N_col, mean = 0,     sd = 0.1), cov = cov(matrix(rnorm(N_col*N_col, 0, 0.5), nrow = N_col, ncol = N_col))),
    cluster3 = list(mean = rnorm(N_col, mean = 0.75,  sd = 0.1), cov = cov(matrix(rnorm(N_col*N_col, 0, 0.5), nrow = N_col, ncol = N_col)))
  ),
  
  # 5. High overlap, full covariance + correlation
  list(
    cluster1 = list(mean = rnorm(N_col, mean = -0.3, sd = 0.1), cov = cov(matrix(rnorm(N_col*N_col, 1, 0.3), nrow = N_col, ncol = N_col))),
    cluster2 = list(mean = rnorm(N_col, mean = 0,    sd = 0.1), cov = cov(matrix(rnorm(N_col*N_col, 1, 0.3), nrow = N_col, ncol = N_col))),
    cluster3 = list(mean = rnorm(N_col, mean = 0.3,  sd = 0.1), cov = cov(matrix(rnorm(N_col*N_col, 1, 0.3), nrow = N_col, ncol = N_col)))
  )
)

param_labels <- c(
  "Well-separated means, identity cov",
  "Closer means, identity cov",
  "Closer still, diag cov varied scale",
  "Full covariance, small overlap",
  "High overlap, full covariance + correlation"
)

# Parallel  
parallel_process <- TRUE

# Generate parameters 
param_set <- params_list[[param_index]] 
param_label <- param_labels[[param_index]]

# Get cores
n_cores <- numCores()

# Start simulation    
data <- simulateGMM(3, n_groups, param_set, n_indiv = 419, n_col = N_col,
                    random_seed = seed,
                    equal_clust = FALSE, equal_groups = FALSE)
true_clusters <- data[[2]]

if (n_groups > 1) {
  true_groups <- data[[3]]
}

data <- data[[1]]

# Classify using GMM
classification <- GMMclassifier(data) # This does the MClust fitting

# Construct empty similarity matrix 
sim_matrix <- matrix(0, nrow = ncol(classification), ncol = ncol(classification),
                     dimnames = list(colnames(classification), colnames(classification)))

# Set seed 
set.seed(seed)

if (parallel_process) {
  
  if (.Platform$OS.type == "windows") {
    cl2 <- makeCluster(n_cores)
    clusterExport(cl2, ls(envir = environment()), envir = environment())
    
    # Parallelised computation of similarity matrix
    sim_mat <- parLapply(cl2, 1:ncol(classification), function(i) {
      sapply(1:ncol(classification), function(j) {
        adjustedRandIndex(classification[[i]], classification[[j]])
      })
    })
    
    # Stop the parallel cluster after the work is done
    stopCluster(cl2)
    
  } else {
    sim_mat <- mclapply(1:ncol(classification), function(i) {
      sapply(1:ncol(classification), function(j) {
        adjustedRandIndex(classification[[i]], classification[[j]])
      })
    }, mc.cores = n_cores)
  }
  
} else {
  sim_mat <- lapply(1:ncol(classification), function(i) {
    sapply(1:ncol(classification), function(j) {
      adjustedRandIndex(classification[[i]], classification[[j]])
    })
  })
}

# Combine the list into a matrix
sim_matrix <- do.call(cbind, sim_mat)

# Convert to dissimilarity matrix
dissim_matrix <- 1 - sim_matrix
dist_mat <- as.dist(dissim_matrix)

# Perform MDS on the distance matrix
mds_result <- cmdscale(dist_mat, k = 2, eig = TRUE)

# Validate MDS output
valid_eigs <- sum(mds_result$eig > 0)
if (valid_eigs < 2 || any(is.na(mds_result$points))) {
  stop("MDS failed: fewer than 2 positive eigenvalues or contains NA values.")
}

# Extract MDS coordinates
mds <- mds_result$points

# Make sure results directory exists
if (!dir.exists("results")) {
  dir.create("results")
}

# Initialize results variable to NULL just in case
results <- NULL

# Fit models to MDS
if (toupper(model) != "OTRIMLE") {
  
  start_time <- Sys.time()
  test <- mclust::mclustBootstrapLRT(mds, modelName = model, nboot = 1000, level = 0.95, maxG = 1)
  end_time <- Sys.time()
  
  if (!is.null(test$p.value)) {
    p <- round(test$p.value, 4)
    results <- data.frame(
      Parameter.ID = param_index,
      Parameter.Label = param_label,
      Columns = N_col,
      `True Number of Groups` = n_groups,
      Method = paste("mclustBootstrapLRT:", model),
      Statistic = paste("p-value:", p),
      `Determined Number of Groups` = ifelse(p < 0.05, ">1", "1"),
      RunTime = paste0(round(as.numeric(end_time - start_time, units = "secs"), 4), "s"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    # Optional: record failure
    results <- data.frame(
      Parameter.ID = param_index,
      Parameter.Label = param_label,
      Columns = N_col,
      `True Number of Groups` = n_groups,
      Method = paste("mclustBootstrapLRT:", model),
      Statistic = "p-value: NA",
      `Determined Number of Groups` = "NA",
      RunTime = NA,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  
  write.csv(results, file = paste0("results/mds_results_", task_id, ".csv"), sep = ",", row.names = FALSE)
  
} else {
  
  results <- tryCatch({
    message("Calling otrimleg()...")
    start_time <- Sys.time()
    model_otrimle <- otrimleg(mds, G = 1:2)
    end_time <- Sys.time()
    
    if (!is.list(model_otrimle)) stop("OTRIMLE did not return a list")
    if (!"ibic" %in% names(model_otrimle)) stop("OTRIMLE result missing 'ibic'")
    
    ibic <- model_otrimle$ibic
    min_ibic <- min(ibic)
    best_G <- which.min(ibic)
    
    data.frame(
      Parameter.ID = param_index,
      Parameter.Label = param_label,
      Columns = N_col,
      `True Number of Groups` = n_groups,
      Method = "OTRIMLE",
      Statistic = paste("min(BIC):", round(min_ibic, 4)),
      `Determined Number of Groups` = as.character(best_G),
      RunTime = paste0(round(as.numeric(end_time - start_time, units = "secs"), 4), "s"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }, error = function(e) {
    message("OTRIMLE error for param set ", param_index, ": ", e$message)
    data.frame(
      Parameter.ID = param_index,
      Parameter.Label = param_label,
      Columns = N_col,
      `True Number of Groups` = n_groups,
      Method = "OTRIMLE",
      Statistic = paste("ERROR:", e$message),
      `Determined Number of Groups` = "NA",
      RunTime = NA,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  
  write.csv(results, file = paste0("results/mds_results_", task_id, ".csv"), row.names = FALSE, sep = ",")
}

