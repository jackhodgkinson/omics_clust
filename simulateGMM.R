## simulateGMM.R
# No external packages needed! 
simulateGMM <- function(n_clust,                                                 # Number of clusters
                        n_groups,                                                # Number of groups of data, 1 by default
                        cluster_params,                                          # List of distribution parameters per cluster
                        n_indiv,                                                 # Number of individuals  
                        n_col,                                                   # Number of columns in simulated data
                        random_seed,                                             # Input random seed for reproducibility
                        cluster_labels = NA,                                     # Input cluster labels, NA by default.
                        equal_clust = TRUE,                                      # If generating cluster labels, ensure each cluster has approx equal individuals
                        equal_groups = TRUE                                      # If n_groups > 1, ensure each group contains approx equal number of cols
                        ){
  
  # Ensure numeric inputs
  n_clust <- as.numeric(n_clust)
  n_groups <- as.numeric(n_groups)
  n_indiv <- as.numeric(n_indiv)
  n_col <- as.numeric(n_col)
  random_seed <- as.numeric(random_seed)

  # Set seed
  set.seed(random_seed)
  
  if (is.na(cluster_labels)){
    # Generate cluster labels for each data point 
    cluster_labs <- seq(1, n_clust)
    
    # Generate random clusters for each n_indiv
    if (equal_clust == FALSE){
    indiv_clust <- sample(cluster_labs, size = n_indiv, replace = TRUE, 
                          prob = {p <- runif(length(cluster_labs)); p / sum(p)})
    }
    else {
      indiv_clust <- sample(cluster_labs, size = n_indiv, replace = TRUE)  
    }
  }
  else {
    indiv_clust <- cluster_labels
  }
  
  # Create empty simulated data matrix
  sim_data <- matrix(NA, ncol = n_col, nrow = n_indiv)
  
  # Generate data 
  for (i in 1:n_col){
    
    # Create vector for the data and params for column i 
    col_data <- numeric(n_indiv)
    mu <- numeric(n_indiv) 
    sigma <- numeric(n_indiv)
    
    # Generate data for each individual based on cluster params
    for (j in 1:n_indiv) {
      
      # Simulate random normals for each protein and cluster
      mu[j] <- cluster_params[[paste0("cluster", indiv_clust[j])]]$mean[i]
      sigma[j] <-  cluster_params[[paste0("cluster", indiv_clust[j])]]$sd[i]
      
      # Generate cluster assignment
      col_data[j] <- rnorm(1, mean = mu[j], sd = sigma[j])
    }
    
    # Add column data to data frame
    sim_data[, i] <- col_data
  }
  
  group <- NULL
  
  # Permute data if different clustering structures required
  if (n_groups > 1) {
    
    # Randomly group 
    if (equal_groups == FALSE){
      repeat {
        p <- runif(n_groups)
        p <- p / sum(p)  # normalize to sum to 1
        group <- sample(1:n_groups, n_col, replace = TRUE, prob = p)
        if (length(unique(group)) == n_groups) break
      }
    }
    
    else {
      group <- sample(1:n_groups, n_col, replace = TRUE)
    }
    
    # Create empty data frame for the clusters for each group 
    group_clusterID <- list()

    # Permute values and clusters within groups 
    for (g in 1:(n_groups - 1)) { 
      permute_cols <- which(group == g)
      order <- sample(n_indiv)
      sim_data[, permute_cols] <- sim_data[order, permute_cols]
      col_name <- paste0("group", g, "_clusterid")
      group_clusterID[[col_name]] <- indiv_clust[order]
    }
  
  # Name the grouping for the final group 
    col_name <- paste0("group", n_groups, "_clusterid")
    group_clusterID[[col_name]] <- indiv_clust
    group_clusterID <- as.data.frame(group_clusterID)
  }

  # Convert to data frame 
  sim_data <- as.data.frame(sim_data)
  
  # Return new dataset 
  if (!is.null(group)){
    return(list("Simulated Data" = sim_data, 
                "Cluster ID per individual per group" = group_clusterID, 
                "Data Group ID" = group))
  }
  else {
    return(list("Simulated Data" = sim_data, 
                "Cluster ID" = indiv_clust))
  }

}