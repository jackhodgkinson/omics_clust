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

# Simulated data 
seed <- 4881
set.seed(seed)
N_col <- 10

# Create empty results dataframe
results <- data.frame(Parameter.ID = numeric(),
                      Columns = numeric(),
                      `True Number of Groups` = numeric(), 
                      Method = character(), 
                      Statistic = character(),
                      `Determined Number of Groups` = character(),
                      RunTime = numeric(),
                      stringsAsFactors = FALSE,
                      check.names = FALSE)

# Parallel 
parallel_process <- TRUE
  
# Specify parameters 
params <- list(
  cluster1 = list(mean = rnorm(N_col, mean = -1.5, sd = 0.1), cov = diag(N_col)),
  cluster2 = list(mean = rnorm(N_col, mean = 0,    sd = 0.1), cov = diag(N_col)),
  cluster3 = list(mean = rnorm(N_col, mean = 1.5,  sd = 0.1), cov = diag(N_col))
)

# Get cores
n_cores <- numCores()

# Simulate data
n_groups <- 2

data <- simulateGMM(3, n_groups, params, n_indiv = 419, n_col = N_col,
                    random_seed = seed,
                    equal_clust = FALSE, equal_groups = FALSE)
true_clusters <- data[[2]]

if (n_groups > 1) {
  true_groups <- data[[3]]
}

data <- data[[1]]

# Plot data
data_plot <- as.matrix(data)
annotationRow <- as.data.frame(true_clusters[1])
names(annotationRow) <- "Clusters"

annotationCol <- as.data.frame(as.factor(true_groups))
names(annotationCol) <- "ProteinClusters"
rownames(annotationCol) <- colnames(data_plot)

rownames(data_plot) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
annotationCol$ProteinClusters <- as.factor(annotationCol$ProteinCluster)
true_clust <- as.numeric(true_clusters[["group1_clusterid"]])
ordered_data <- data_plot[order(true_clust), order(annotationCol$ProteinClusters)]
pheatmap::pheatmap(
  ordered_data,
  annotation_row = annotationRow,
  annotation_col = annotationCol[order(annotationCol$ProteinClusters), , drop = FALSE],
  cluster_rows = FALSE
)

# Get distance matrix
classification <- GMMclassifier(data) # This does the MClust fitting - can replace with LCMM

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

# Pheatmap
pheatmap::pheatmap(sim_matrix)

# Convert to dissimilarity matrix
dissim_matrix <- (1 - sim_matrix) + abs(matrix(rnorm(N_col*N_col, 0, 0.001), N_col , N_col))
# add noise for easier separation in MDS subspace.

#dissim_matrix <- dissim_matrix2 
dist_mat <- as.dist(dissim_matrix)

# Perform MDS on the distance matrix
mds_result <- cmdscale(dist_mat, k = 2, eig = TRUE)

# Validate MDS output
valid_eigs <- sum(mds_result$eig > 0)
if (valid_eigs < 2 || any(is.na(mds_result$points))) {
  stop("MDS failed: fewer than 2 positive eigenvalues or contains NA values.")
}

# Extract mds points
mds <- as.data.frame(mds_result$points)

# Plot MDS output
ggplot(mds, aes(x = mds[,1], y = mds[,2])) + 
  geom_point() + 
  labs(x = "MDS1",
       y = "MDS2",
       title = "2D Multi-dimensional scaling plot of dissimilarity matrix")

# Fit Gaussian and Gaussian Mixture (2 components) to MDS
for (m in models) {
  start_time <- Sys.time()
  test <- mclust::mclustBootstrapLRT(mds, modelName = "VII", nboot = 1000, level = 0.95, maxG = 1)
  if (test$p.value < 0.05) {
    g = 2:N_col 
  } else {
    g = 1
    }
  test2 <- mclust::Mclust(mds, G = g)
  plot(test2, "classification")
  end_time <- Sys.time()
  if (!is.null(test$p.value)) {
    p <- round(test$p.value, 4)
    results <- rbind(results, data.frame(
      #Parameter.ID = param_index,
      #Parameter.Label = param_label,
      Columns = N_col,
      `True Number of Groups` = n_groups,
      Method = paste("mclustBootstrapLRT:", m),
      Statistic = paste("p-value:", p),
      `Determined Number of Groups` = ifelse(p < 0.05, ">1", "1"),
      RunTime = end_time - start_time
    ))
  } else {
    # Optional: record failure
    results <- rbind(results, data.frame(
      #Parameter.ID = param_index,
      #Parameter.Label = param_label,
      Columns = N_col,
      `True Number of Groups` = n_groups,
      Method = paste("mclustBootstrapLRT:", m),
      Statistic = "p-value: NA",
      `Determined Number of Groups` = "NA",
      RunTime = "NA"
    ))
  }
#}

# OTRIMLE
otrimle_result <- tryCatch({
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
    #Parameter.ID = param_index,
    #Parameter.Label = param_label,
    Columns = N_col,
    `True Number of Groups` = n_groups,
    Method = "OTRIMLE",
    Statistic = paste("min(BIC):", toString(round(min_ibic, 4))),
    `Determined Number of Groups` = as.character(which.min(ibic)),
    RunTime = Sys.time() - start_time
  )
  
  g <- which.min(ibic)
  plot_otrimle <- otrimle(mds, G=g)
  plot(plot_otrimle, "clustering", mds)
  
}, error = function(e) {
  message("OTRIMLE error : ", e$message)
  data.frame(
    #Parameter.ID = param_index,
    #Parameter.Label = param_label,
    Columns = N_col,
    `True Number of Groups` = n_groups,
    Method = "OTRIMLE",
    Statistic = paste("ERROR:", e$message),
    `Determined Number of Groups` = "NA",
    RunTime = NA
  )
})

results <- rbind(results, otrimle_result)


#write.csv(results, file = "mds_results.csv", row.names = FALSE)