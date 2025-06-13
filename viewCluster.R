# viewCluster.R 
viewCluster <- function(classification,
                        parallel = TRUE 
                        ) {
  
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
  dissim_matrix <- 1 - sim_matrix
  dist_mat <- as.dist(dissim_matrix)
  
  # Perform MDS on the distance matrix
  mds_result <- cmdscale(dist_mat, k = 2, eig = TRUE)
  mds <- mds_result$points
  mds_plot_data <- as.data.frame(mds)
  
  mds_plot <- ggplot(mds_plot_data, aes(x = V1, y = V2)) + 
    geom_point() + 
    labs(
      x = "MDS1",
      y = "MDS2",
      title = "Multi-dimensional Scaling (MDS) plot of the ARI dissimilarity matrix"
    ) 
  
  ## Hierarchcial Clustering on views
  # Hierarchical Clustering 
  hclust <- hclust(dist_mat, method = linkage)
  opt <- NbClust(diss = dist_mat, distance = NULL, method = linkage, 
                 min.nc = 2, max.nc = floor(sqrt(ncol(classification))), index = index)
  
  ## Split data into different views from dendrogram
  # Choose optimal number of clusters
  if linkage %in% c("mcclain","cindex") {
    if (is.vector(opt$Best.nc)) {
      opt_ng <- as.numeric(opt$Best.nc[1])
    } else {
      opt_ng <- as.numeric(names(which.min(table(opt$Best.nc[1, ]))))
    }
  }
  
  else {
    if (is.vector(opt$Best.nc)) {
      opt_ng <- as.numeric(opt$Best.nc[1])
    } else {
      opt_ng <- as.numeric(names(which.max(table(opt$Best.nc[1, ]))))
    }
  }
  
}