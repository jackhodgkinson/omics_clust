## simulateGMM.R
# MASS needed
simulateGMM <- function(n_clust,                                                 # Number of clusters
                        n_groups,                                                # Number of groups of data, 1 by default
                        cluster_params,                                          # List of distribution parameters per cluster
                        n_indiv,                                                 # Number of individuals  
                        n_col,                                                   # Number of columns in simulated data
                        random_seed,                                             # Input random seed for reproducibility
                        cluster_labels = NULL,                                   # Input cluster labels, NULL by default.
                        group_labels = NULL,                                     # Input group labels, NA by default
                        equal_clust = TRUE,                                      # If generating cluster labels, ensure each cluster has approx equal individuals
                        equal_groups = TRUE,                                     # If n_groups > 1, ensure each group contains approx equal number of cols
                        _parallel = FALSE                                          # Use parallel processsing. Default is TRUE
                        ){
  
  # Load relevant packages
  library(MASS)
  library(parallel)
  
  # Ensure numeric inputs
  n_clust <- as.numeric(n_clust)
  n_groups <- as.numeric(n_groups)
  n_indiv <- as.numeric(n_indiv)
  n_col <- as.numeric(n_col)
  random_seed <- as.numeric(random_seed)

  # Set seed
  set.seed(random_seed)
  
  # Validate parameters
  if (length(cluster_params) != n_clust) {
    stop(paste0("The number of clusters (", n_clust, ") does not match the number of elements in cluster_params (", length(cluster_params), ")."))
  }
  invisible(lapply(seq_len(n_clust), function(k) {
    clust_name <- paste0("cluster", k)
    params <- cluster_params[[clust_name]]
    if (is.null(params)) stop(paste("Missing parameters for", clust_name))
    if (is.null(params$mean)) stop(paste("Missing 'mean' for", clust_name))
    if (length(params$mean) != n_col) stop(paste("Length of 'mean' must match n_col for", clust_name))
    if (is.null(params$cov) && is.null(params$sd)) {
      stop(paste("Must provide either 'cov' or 'sd' for", clust_name))
    }
    if (!is.null(params$sd) && length(params$sd) != n_col) {
      stop(paste("Length of 'sd' must match n_col for", clust_name))
    }
    if (!is.null(params$cov) && (!is.matrix(params$cov) || any(dim(params$cov) != n_col))) {
      stop(paste("Covariance matrix for", clust_name, "must be", n_col, "x", n_col))
    }
    if (!is.null(params$sd) && !is.null(params$cov)) {
      warning(paste("Both 'sd' and 'cov' provided for", clust_name, "- using 'cov'"))
    }
  }))
  
  # Generate cluster labels
  indiv_clust <- if (is.null(cluster_labels)){
    if (!equal_clust){
      sample(seq_len(n_clust), size = n_indiv, replace = TRUE, 
             prob = {p <- runif(n_clust); p / sum(p)})
    } else {
      sample(seq_len(n_clust), size = n_indiv, replace = TRUE)  
    } 
  } else cluster_labels
  
  # Data simulation function
  sim_clust_data <- function (k) {
    idx <- which(indiv_clust == k)
    if (length(idx) == 0) return(NULL)
    param <- cluster_params[[paste0("cluster",k)]]
    mu <- param$mean
    Sigma <- if (!is.null(param$cov)) param$cov else diag(param$sd^2)
    data <- mvrnorm(n = length(idx), mu = mu, Sigma = Sigma)
    list(index = idx, data = data)
  }
  
  # Generate data 
  if (_parallel && n_indiv > 1000 && n_clust > 3) {
    if (.Platform$OS.type != "windows"){
      cluster_results <- mclapply(seq_len(n_clust), sim_clust_data, 
                                  mc.cores = detectCores() - 1)
    }
    else {
      cl <- makeCluster(detectCores() - 1)
      clusterExport(cl, ls(envir = environment()), envir = environment())
      cluster_results <- parLapply(cl, seq_len(n_clust), sim_clust_data)
      stopCluster(cl)
    }
  } else { 
      cluster_results <- lapply(seq_len(n_clust), sim_clust_data)
    }
      
    # Create data matrix 
    sim_data <- matrix(NA, nrow = n_indiv, ncol = n_col)
    invisible(lapply(cluster_results, function(res) {
      if (!is.null(res)) sim_data[res$index, ] <<- res$data
    }))
    
    # Assign empty group
    group <- NULL
    
    # Permute data if different clustering structures required
    if (n_groups > 1) {
      group <- if (!is.null(group_labels)) {
        group_labels
      } else if (equal_groups == FALSE & !is.null(group)){
        repeat {
          p <- runif(n_groups)
          p <- p / sum(p)  # normalize to sum to 1
          group <- sample(seq_len(n_groups), n_col, replace = TRUE, prob = p)
          if (length(unique(group)) == n_groups) break
        }
        group 
      } else {
        sample(seq_len(n_groups), n_col, replace = TRUE)
      }
      
      # Run in parallel if needed
      if (_parallel && n_indiv > 1000 && n_groups > 2) {
        if (.Platform$OS.type != "windows") {
          cl <- mclapply(seq_len(n_groups - 1), function(g) {
            permute_cols <- which(group == g)
            order <- sample(n_indiv)
            sim_data[, permute_cols] <<- sim_data[order, permute_cols]
            out <- indiv_clust[order]
            names(out) <- paste0("group", g, "_clusterid")
            out
          }, mc.cores = detectCores() - 1)
        } else {
          cl <- makeCluster(n_cores)
          clusterExport(cl, c("sim_data", "group", "n_indiv", "indiv_clust"), envir = environment())
          permuted_list <- parLapply(cl, seq_len(n_groups - 1), function(g) {
            permute_cols <- which(group == g)
            order <- sample(n_indiv)
            sim_data[, permute_cols] <<- sim_data[order, permute_cols]
            out <- indiv_clust[order]
            names(out) <- paste0("group", g, "_clusterid")
            out
          })
          stopCluster(cl)
        }
    } else {
      # Serial column permutations for groups
      permuted_list <- lapply(seq_len(n_groups - 1), function(g) {
        permute_cols <- which(group == g)
        order <- sample(n_indiv)
        sim_data[, permute_cols] <<- sim_data[order, permute_cols]
        out <- indiv_clust[order]
        names(out) <- paste0("group", g, "_clusterid")
        out
      })
    }
      
    names(permuted_list) <- paste0("group", seq_len(n_groups-1), "_clusterid")
    permuted_list[[paste0("group", n_groups, "_clusterid")]] <- indiv_clust
    group_clusterID <- as.data.frame(permuted_list)
    
    sim_data <- as.data.frame(sim_data)
    
    return(list("Simulated Data" = sim_data,
                "Cluster ID per individual per group" = group_clusterID,
                "Group ID"= group))
  } else {
      sim_data <- as.data.frame(sim_data)
      return(list("Simulated Data" = sim_data,
                  "Cluster ID" = indiv_clust))
  }
}
