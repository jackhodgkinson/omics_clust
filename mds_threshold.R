## MDS Thresholding 
# Turn off graphics
#graphics.off()

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
N_col <- c(10, 100, 500, 1000)

# Create empty results dataframe
results <- data.frame(Columns = numeric(),
                      `True Number of Groups` = numeric(), 
                      Method = character(), 
                      Statistic = character(),
                      `Determined Number of Groups` = character(),
                      RunTime = numeric(),
                      stringsAsFactors = FALSE,
                      check.names = FALSE)

# Parallel 
parallel_process <- TRUE

# Loop over different values of proteins
for (i in N_col) {
  
  # Less defined clusters with correlation
  params <- list(
    cluster1 = list(mean = rnorm(i, mean = -1, sd = 0.1), cov = cov(matrix(rnorm(i*i, mean = 0, sd = 0.5), nrow = i, ncol = i))),
    cluster2 = list(mean = rnorm(i, mean = 0, sd = 0.1), cov = cov(matrix(rnorm(i*i, mean = 1, sd = 0.25), nrow = i, ncol = i))),
    cluster3 = list(mean = rnorm(i, mean = 1, sd = 0.1), cov = cov(matrix(rnorm(i*i, mean = -1, sd = 0.15), nrow = i, ncol = i)))
  )
  
  # Simulate data
  n_groups <- c(1, 2, 3, 4)
  
  # Get cores
  n_cores <- numCores()
  
  for (n in n_groups) {
    
    data <- simulateGMM(3, n, params, n_indiv = 419, n_col = i,
                         random_seed = seed,
                         equal_clust = FALSE, equal_groups = FALSE)
    true_clusters <- data[[2]]
    
    if (n > 1) {
      true_groups <- data[[3]]
    }
    
    data <- data[[1]]
    
    # Plot data
    # data_plot <- as.matrix(data)
    # annotationRow <- as.data.frame(true_clusters[1])
    # names(annotationRow) <- "Clusters"
    # 
    # annotationCol <- as.data.frame(as.factor(true_groups))
    # names(annotationCol) <- "ProteinClusters"
    # rownames(annotationCol) <- colnames(data_plot)
    # 
    # rownames(data_plot) <- rownames(annotationRow)
    # annotationRow$Clusters <- as.factor(annotationRow$Clusters)
    # annotationCol$ProteinClusters <- as.factor(annotationCol$ProteinCluster)
    # true_clust <- as.numeric(true_clusters[["group1_clusterid"]])
    # ordered_data <- data_plot[order(true_clust), order(annotationCol$ProteinClusters)]
    # pheatmap::pheatmap(
    #   ordered_data,
    #   annotation_row = annotationRow,
    #   annotation_col = annotationCol[order(annotationCol$ProteinClusters), , drop = FALSE],
    #   cluster_rows = FALSE
    # )
    
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
    }
    
    else {
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
    mds <- cmdscale(dist_mat)
    
    # Fit Gaussian and Gaussian Mixture (2 components) to MDS
    models <- c("EII","VII","EEI","VEI","EVI","VVI","EEE",
                "VEE","EVE","VVE","EEV","VEV","EVV","VVV")
    for (m in models) {
      start_time <- Sys.time()
      test <- mclust::mclustBootstrapLRT(mds, modelName = m, nboot = 1000, level = 0.95, maxG = 1)
      end_time <- Sys.time()
      if (!is.null(test$p.value)) {
        p <- round(test$p.value, 4)
        results <- rbind(results, data.frame(
          Columns = i,
          `True Number of Groups` = n,
          Method = paste("mclustBootstrapLRT:", m),
          Statistic = paste("p-value:", p),
          `Determined Number of Groups` = ifelse(p < 0.05, ">1", "1"),
          RunTime = end_time - start_time
        ))
      } else {
        # Optional: record failure
        results <- rbind(results, data.frame(
          Columns = i,
          `True Number of Groups` = n,
          Method = paste("mclustBootstrapLRT:", m),
          Statistic = "p-value: NA",
          `Determined Number of Groups` = "NA",
          RunTime = "NA"
        ))
      }
    }
    
    # OTRIMLE
    start_time = Sys.time()
    model_otrimle <- otrimleg(mds, G=1:2)
    end_time = Sys.time()
    bic <- min(model_otrimle$ibic)
    
    results <- rbind(results, data.frame(Columns = i,`True Number of Groups` = n,
                                         Method = "OTRIMLE",
                                         Statistic = paste("min(BIC):", bic),
                                         `Determined Number of Groups` = ifelse(which.min(model_otrimle$ibic) == 1, "1", ">1"),
                                         RunTime = end_time - start_time)
    )
  }
}

write.csv(results, file = "mds_results.csv", row.names = FALSE)

