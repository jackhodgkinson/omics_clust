# simulateLCMM
simulateLCMM <- function(subject_data = NULL,                                     # Provide subject data. Default is NULL and function will generate subject data.
                         timepoints,                                              # Timepoints for longitudinal data
                         n_clust,                                                 # Number of clusters
                         n_groups,                                                # Number of groups of data, 1 by default
                         cluster_params,                                          # List of distribution parameters per cluster
                         n_indiv,                                                 # Number of individuals
                         n_col,                                                   # Number of columns in simulated data
                         random_seed,                                             # Input random seed for reproducibility
                         missing = FALSE,                                         # Boolean. Default is FALSE.
                         missing_perc = 0.05,                                     # Percentage of missing data if missing == TRUE. Default is 0.05.
                         missing_timepoints = NULL,                               # The index of timepoints you would like to add missing data to. Default is NULL, so selection is random.
                         timepoint_perc = 0.1,                                    # The proportion of participants missing 1:n-1 timepoints. Can be a scalar or a vector of length n-1 timepoints. Default is 0.1
                         timepoint_noise = TRUE,                                  # Add noise to the timepoints. Default is TRUE.
                         timepoint_sd = 0,                                        # Specify the standard deviation of the timepoints to add noise. Can be a single number of a vector the same length as timepoints. Default is 0.
                         cluster_labels = NULL,                                   # Input cluster labels, NULL by default.
                         group_labels = NULL,                                     # Input group labels, NULL by default
                         equal_clust = TRUE,                                      # Boolean. If generating cluster labels, ensure each cluster has approx equal individuals. Default is TRUE.
                         equal_groups = TRUE,                                     # Boolean. If n_groups > 1, ensure each group contains approx equal number of columns. Default is TRUE.
                         parallel_process = TRUE                                  # Boolean. Use parallel processsing. Default is TRUE
                         ){
  
  # ==== PACAKGES ==== 
  library(MASS)
  library(parallel)
  source("numCores.R")

  # ==== PARAMETER VALIDATION ====
  if (!is.numeric(timepoints) || length(timepoints) < 1) stop("`timepoints` must be a numeric vector.")
  if (!is.numeric(n_clust) || n_clust < 1) stop("`n_clust` must be a positive integer.")
  if (!is.numeric(n_groups) || n_groups < 1) stop("`n_groups` must be a positive integer.")
  if (!is.list(cluster_params)) stop("`cluster_params` must be a list.")
  if (!is.numeric(n_indiv) || n_indiv < 1) stop("`n_indiv` must be a positive integer.")
  if (!is.numeric(n_col) || n_col < 1) stop("`n_col` must be a positive integer.")
  if (!is.numeric(random_seed)) stop("`random_seed` must be numeric.")
  if (!is.logical(missing)) stop("`missing` must be TRUE or FALSE.")
  
  if (!is.null(missing_timepoints)) {
    if (!all(missing_timepoints %in% timepoints)) {
      stop("All `missing_timepoints` must be in `timepoints`.")
    }
  } else {
    if (!is.numeric(timepoint_perc) || any(timepoint_perc < 0)) {
      stop("timepoint_perc must be numeric and non-negative")
    }
  }
  
  if (length(missing_perc) != 1 && !is.null(missing_timepoints) &&
      length(missing_perc) != length(missing_timepoints)) {
    stop("`missing_perc` must be a single value or match length of `missing_timepoints`.")
  }
  
  if (any(missing_perc < 0 | missing_perc > 1)) {
    stop("`missing_perc` values must be between 0 and 1.")
  }
  
  # ==== GENERATE SUBJECT DATA ====  
  if (is.null(subject_data)) {     
    data_hlme <- data.frame(
      Subject_ID = rep(1:n_indiv, each = length(timepoints)),
      Time = rep(timepoints, times = n_indiv)
    )   
  } else {     
    required_cols <- c("Subject_ID", "Time")     
    missing_cols <- setdiff(required_cols, colnames(subject_data))     
    if (length(missing_cols) > 0) {       
      stop(paste("Subject data is missing required column(s):", paste(missing_cols, collapse = ", ")))     
    }     
    expected_rows <- n_indiv * length(timepoints)     
    if (nrow(subject_data) != expected_rows) {       
      warning(paste("Subject data has", nrow(subject_data), "rows; expected", expected_rows))     
    }     
    data_hlme <- subject_data   
  }

  # ==== APPLY MISSINGNESS ====
  if (missing) {
    set.seed(random_seed) 
    subj <- unique(data_hlme$Subject_ID)
    
    n_miss <- numeric(length(timepoint_perc))
    
    # Calculate number of participants missing at least n-1,...,1 timepoint
    for (k in seq_along(timepoint_perc)) {
        n_miss[k] <- round(timepoint_perc[k] * length(subj))
    }
    n_miss <- c("0" = length(subj) - sum(n_miss), setNames(n_miss, seq_along(n_miss)))
    
    # Initialise the target of missing counts per timepoint
    target_miss_tp <- round(missing_perc * length(subj))
    names(target_miss_tp) <- timepoints
    
    # Extend n_miss to the number of participants and permute randomly. 
    n_miss <- rep(as.integer(names(n_miss)), times = n_miss)
    set.seed(random_seed) 
    n_miss <- sample(n_miss, length(n_miss), replace = FALSE)
    
    # Matrix
    mat <- matrix(runif(length(n_miss * )))
    
    # Matrix of missingness
    missing_mat <- matrix(FALSE, nrow = length(subj), ncol = length(timepoints))
    
    
    # Loop over each timepoint and mark n_miss timepoints as missing using target_miss_tp
    for (i in seq_along(subj)){
      if (n_miss[i] == 0) next
      tps <- names(which(colSums(missing_mat) < target_miss_tp))
      
      n_assign <- min(n_miss[i], length(tps))
      
      tp_assign <- sample(tps, n_assign)
      
      missing_mat[i, tp_assign] <- TRUE
    }
    
    return(missing_mat)
  }
  
  # if (missing) {
  #   set.seed(random_seed)
  #   all_subjects <- unique(data_hlme$Subject_ID)
  #   total_participants <- length(all_subjects)
  #   
  #   # Define hierarchical missingness proportions (must sum <= 1)
  #   prop_missing_3 <- 0.0098
  #   prop_missing_2 <- 0.0802
  #   prop_missing_1 <- 0.3521
  #   
  #   # Calculate number of participants missing at least 3, 2, and 1 timepoints
  #   n_miss_3 <- round(prop_missing_3 * total_participants)
  #   n_miss_2 <- round(prop_missing_2 * total_participants)
  #   n_miss_1 <- round(prop_missing_1 * total_participants)
  #   
  #   # Sample participants hierarchically (no overlap between groups)
  #   miss_3_participants <- sample(all_subjects, n_miss_3, replace = FALSE)
  #   remaining_for_2 <- setdiff(all_subjects, miss_3_participants)
  #   miss_2_participants <- sample(remaining_for_2, n_miss_2, replace = FALSE)
  #   remaining_for_1 <- setdiff(all_subjects, c(miss_3_participants, miss_2_participants))
  #   miss_1_participants <- sample(remaining_for_1, n_miss_1, replace = FALSE)
  #   
  #   # Initialize missing participants per timepoint list
  #   missing_participants_per_tp <- vector("list", length(timepoints))
  #   names(missing_participants_per_tp) <- as.character(timepoints)
  #   for (tp in timepoints) {
  #     missing_participants_per_tp[[as.character(tp)]] <- integer(0)
  #   }
  #   
  #   # Assign missing timepoints to 3-missing group
  #   for (subj in miss_3_participants) {
  #     tps <- sample(timepoints, 3, replace = FALSE)
  #     for (tp in tps) {
  #       missing_participants_per_tp[[as.character(tp)]] <- c(missing_participants_per_tp[[as.character(tp)]], subj)
  #     }
  #   }
  #   
  #   # Assign missing timepoints to 2-missing group
  #   for (subj in miss_2_participants) {
  #     tps <- sample(timepoints, 2, replace = FALSE)
  #     for (tp in tps) {
  #       missing_participants_per_tp[[as.character(tp)]] <- c(missing_participants_per_tp[[as.character(tp)]], subj)
  #     }
  #   }
  #   
  #   # Assign missing timepoints to 1-missing group
  #   for (subj in miss_1_participants) {
  #     tp <- sample(timepoints, 1)
  #     missing_participants_per_tp[[as.character(tp)]] <- c(missing_participants_per_tp[[as.character(tp)]], subj)
  #   }
  #   
  #   # Remove duplicates inside each timepoint list
  #   for (tp in timepoints) {
  #     missing_participants_per_tp[[as.character(tp)]] <- unique(missing_participants_per_tp[[as.character(tp)]])
  #   }
  #   
  #   # Collect rows to remove based on missing participants per timepoint
  #   rows_to_remove <- integer(0)
  #   for (tp in timepoints) {
  #     tp_char <- as.character(tp)
  #     subs <- missing_participants_per_tp[[tp_char]]
  #     for (subj in subs) {
  #       rows <- which(data_hlme$Subject_ID == subj & data_hlme$Time == tp)
  #       rows_to_remove <- c(rows_to_remove, rows)
  #     }
  #   }
  #   
  #   # Remove rows to create missingness
  #   data_hlme <- data_hlme[-rows_to_remove, ]
  #   
  #   # --- Check actual missing percentages ---
  #   missing_counts <- sapply(timepoints, function(tp) {
  #     length(unique(setdiff(all_subjects, unique(data_hlme$Subject_ID[data_hlme$Time == tp]))))
  #   })
  #   actual_missing_perc <- missing_counts / total_participants
  #   
  #   # Print summary
  #   cat("Target missing percentages per timepoint:\n")
  #   print(setNames(round(c(prop_missing_1, prop_missing_2, prop_missing_3), 4), c("1_tp", "2_tp", "3_tp")))
  #   
  #   cat("\nActual missing participants per timepoint:\n")
  #   print(setNames(missing_counts, as.character(timepoints)))
  #   
  #   cat("\nActual missing percentages per timepoint:\n")
  #   print(round(actual_missing_perc, 4))
  #   
  #   cat("\nDifference from target missing_perc (input vector):\n")
  #   print(round(actual_missing_perc - missing_perc, 4))
  # }
  # 
  # ==== ADD TIMEPOINT NOISE ====
  non_missing_idx <- which(!is.na(data_hlme$Time))

  if (timepoint_noise) {
    if (length(timepoint_sd) == 1) {
      sd_vals <- rep(timepoint_sd, length(non_missing_idx))
    } else if (length(timepoint_sd) == length(non_missing_idx)) {
      sd_vals <- timepoint_sd
    } else {
      stop("`timepoint_sd` must be a scalar or a vector matching the number of non-missing timepoints.")
    }

    data_hlme$Time[non_missing_idx] <- data_hlme$Time[non_missing_idx] +
      rnorm(length(non_missing_idx), mean = 0, sd = sd_vals)
  }

  # ==== CLUSTER LABELS ====
  set.seed(random_seed)
  indiv_clust <- if (is.null(cluster_labels)) {
    if (!equal_clust) {
      sample(seq_len(n_clust), size = n_indiv, replace = TRUE,
             prob = {p <- runif(n_clust); p / sum(p)})
    } else {
      rep_len(rep(1:n_clust, length.out = n_indiv), n_indiv)
    }
  } else cluster_labels

  indiv_clust_long <- rep(NA_integer_, nrow(data_hlme))
  subject_ids <- data_hlme$Subject_ID
  for (i in seq_along(indiv_clust)) {
    subj_id <- i  # Assuming Subject_ID = 1:n_indiv
    rows_subj <- which(subject_ids == subj_id)
    indiv_clust_long[rows_subj] <- indiv_clust[i]
  }
  #

  # ==== SIMULATION FUNCTION ====

  n_cores <- numCores()
  sim_longitud_data <- function(params, subject_data, ID, indiv_clust, n_col, parallel_proc = TRUE) {
    # ==== Parameter validation ====
    if (!is.list(params) || length(params) == 0) {
      stop("params must be a non-empty list of cluster parameter lists.")
    }

    for (k in seq_along(params)) {
      cluster_name <- paste0("cluster", k)
      if (!cluster_name %in% names(params)) {
        stop(paste("params missing element:", cluster_name))
      }

      p <- params[[cluster_name]]

      if (!("fixed_intercept" %in% names(p)) || length(p$fixed_intercept) != n_col) {
        stop(paste(cluster_name, ": fixed_intercept must be vector of length", n_col))
      }
      if (!("fixed_slope" %in% names(p)) || length(p$fixed_slope) != n_col) {
        stop(paste(cluster_name, ": fixed_slope must be vector of length", n_col))
      }

      if (!("resid_sd" %in% names(p)) || length(p$resid_sd) != n_col) {
        stop(paste(cluster_name, ": resid_sd must be vector of length", n_col))
      }
      if (any(p$resid_sd < 0)) {
        stop(paste(cluster_name, ": resid_sd values must be non-negative"))
      }

      if (!("random_cov" %in% names(p)) || !is.matrix(p$random_cov)) {
        stop(paste(cluster_name, ": random_cov must be a matrix"))
      }

      n_re <- nrow(p$random_cov)
      if (!all(dim(p$random_cov) == c(n_re, n_re))) {
        stop(paste(cluster_name, ": random_cov must be a square matrix"))
      }
      if (!isSymmetric(p$random_cov)) {
        stop(paste(cluster_name, ": random_cov must be symmetric"))
      }
      eigenvalues <- eigen(p$random_cov, symmetric = TRUE)$values
      if (any(eigenvalues <= 0)) {
        stop(paste(cluster_name, ": random_cov must be positive definite"))
      }
    }

    # Ensure indiv_clust has names matching subject IDs
    if (is.null(names(indiv_clust))) {
      unique_ids <- unique(subject_data[[ID]])
      if(length(unique_ids) == length(indiv_clust)) {
        names(indiv_clust) <- as.character(unique_ids)
      } else {
        stop("Length of indiv_clust does not match unique subject IDs in subject_data, and no names found.")
      }
    }

    sim_cluster <- function(k) {
      param <- params[[paste0("cluster", k)]]
      n_re <- nrow(param$random_cov)

      cluster_subjects <- names(indiv_clust)[indiv_clust == k]
      if (length(cluster_subjects) == 0) return(NULL)

      sim_mat <- matrix(NA, nrow = nrow(subject_data), ncol = n_col)

      for (subj_id in cluster_subjects) {
        subject_rows <- which(subject_data[[ID]] == subj_id & !is.na(subject_data$Time))
        if (length(subject_rows) == 0) next

        time_vals <- subject_data$Time[subject_rows]

        b_i <- MASS::mvrnorm(1, mu = rep(0, n_re), Sigma = param$random_cov)

        outcome_mat <- sapply(seq_len(n_col), function(col_idx) {
          fixed_int <- param$fixed_intercept[col_idx]
          fixed_slo <- param$fixed_slope[col_idx]
          resid_sd <- param$resid_sd[col_idx]

          # Generate separate random effects per feature
          b_i_col <- c(rnorm(1, 0, sqrt(param$random_cov[1,1])),
                                  if (n_re >= 2) rnorm(1, 0, sqrt(param$random_cov[2,2])) else NULL)
          re_contrib <- b_i_col[1] + if (n_re >= 2) b_i_col[2] * time_vals else 0

          mu_vec <- fixed_int + fixed_slo * time_vals + re_contrib
          rnorm(length(time_vals), mean = mu_vec, sd = resid_sd)
        })

        sim_mat[subject_rows, ] <- outcome_mat
      }

      sim_mat
    }

    if (parallel_proc) {
      # Detect OS for parallel backend
      os_type <- .Platform$OS.type
      clus <- NULL

      if (os_type == "windows") {
        clus <- parallel::makeCluster(n_cores)
        # Export needed variables and libraries
        parallel::clusterExport(clus, varlist = c("params", "subject_data", "ID", "indiv_clust", "n_col", "sim_cluster"), envir = environment())
        parallel::clusterEvalQ(clus, library(MASS))
        sim_results <- parallel::parLapply(clus, seq_along(params), sim_cluster)
        parallel::stopCluster(clus)
      } else {
        # Unix/macOS: use mclapply
        sim_results <- parallel::mclapply(seq_along(params), sim_cluster, mc.cores = n_cores)
      }
    } else {
      sim_results <- lapply(seq_along(params), sim_cluster)
    }

    sim_data <- matrix(NA, nrow = nrow(subject_data), ncol = n_col)
    for (k in seq_along(sim_results)) {
      res <- sim_results[[k]]
      if (!is.null(res)) {
        non_na_rows <- which(rowSums(!is.na(res)) > 0)
        sim_data[non_na_rows, ] <- res[non_na_rows, ]
      }
    }

    colnames(sim_data) <- paste0("V", seq_len(n_col))
    result <- cbind(subject_data, sim_data)
    return(result)
  }

  # ==== RUN SIMULATION ====
  sim_data <- sim_longitud_data(cluster_params, data_hlme, "Subject_ID",
                               indiv_clust_long, n_col, parallel_proc = parallel_process)

  # ==== GENERATE VIEWS ====
  # Initialise empty group vector
  group <- NULL

  # Generate group labels if n_groups>1
  if (n_groups > 1) {
    if (!is.null(group_labels)) {
      group <- group_labels
    } else if (equal_groups == FALSE & !is.null(group)) {
      repeat {
        p <- runif(n_groups)
        p <- p / sum(p)
        group <- sample(seq_len(n_groups), N_col, replace = TRUE, prob = p)
        if (length(unique(group)) == n_groups) break
      }
    } else {
      group <- sample(seq_len(n_groups), N_col, replace = TRUE)
    }

    # Split the dataset into the different groups
    group_list <- list()
    split_data <- sim_data %>%
      dplyr::select(-Subject_ID, -Time)

    unique_groups <- unique(group)
    for (g in unique_groups[-length(unique_groups)]) {  # Exclude last group
      view <- split_data[, group == g, drop = FALSE]
      view <- cbind(sim_data[, c("Subject_ID", "Time")], view)
      name <- paste0("Group", g)
      group_list[[name]] <- view
    }

    # Empty list for the permuted groups
    group_permute <- list()

    # Permute the columns per group
    for (i in seq_along(names(group_list))) {
      name <- names(group_list)[i]
      group_df <- group_list[[name]]

      permute_col <- function(data, ID, cluster, random_seed) {
        set.seed(random_seed)
        clust_data <- as.data.frame(cbind(cluster, unique(data$Subject_ID)))
        colnames(clust_data)[2] <- "Subject_ID"
        data_clust <- data %>%
          left_join(clust_data, by = "Subject_ID")

        subject_counts <- data_clust %>%
          group_by(Subject_ID) %>%
          summarise(n_timepoints = n())

        data_clust_counts <- data_clust %>%
          left_join(subject_counts, by = "Subject_ID")

        data_by_timepoint <- split(data_clust_counts, data_clust_counts$n_timepoints)
        names(data_by_timepoint) <- paste(seq_along(data_by_timepoint),"Timepoint(s)")
        data_by_timepoint <- lapply(data_by_timepoint, function(df) {df$n_timepoints <- NULL
        return(df)
        })

        data_by_timepoint <- lapply(data_by_timepoint, function(df) {

          orig_ID <- df$Subject_ID
          orig_time <- df$Time
          permute_ID <- sample(unique(df$Subject_ID))
          df <- df %>%
            arrange(match(Subject_ID, permute_ID), Time) %>%
            mutate(Subject_ID = orig_ID,
                   Time = orig_time)
          df
        })

        data2 <- bind_rows(data_by_timepoint) %>%
          arrange(Subject_ID, Time)

        return(data2)
      }

      seed <- random_seed + i
      group_permute[[name]] <- permute_col(group_df, "Subject_ID", indiv_clust, seed)
      names(group_permute[[name]])[names(group_permute[[name]]) == "cluster"] <- paste0(name, "_clusterID")
    }

    # Obtain data for non-permuted group
    non_permute_group_num <- unique_groups[length(unique_groups)]
    non_permute_view <- split_data[, group == non_permute_group_num, drop = FALSE]
    non_permute_group <- cbind(sim_data[, c("Subject_ID", "Time")], non_permute_view)
    non_perm_clust_data <- data.frame(
      Subject_ID = unique(sim_data$Subject_ID),
      cluster = indiv_clust
    )
    non_perm_data_clust <- non_permute_group %>%
      left_join(non_perm_clust_data, by = "Subject_ID")
    names(non_perm_data_clust)[names(non_perm_data_clust) == "cluster"] <- paste0("Group",non_permute_group_num, "_clusterID")

    # List of dataframes
    name <- paste0("Group", non_permute_group_num)
    named_non_perm <- setNames(list(non_perm_data_clust), name)
    dfs <- c(named_non_perm, group_permute)

    # Join all dataframes avoiding duplication of Subject_ID and Time
    final_df <- reduce(dfs, full_join, by = intersect(names(dfs[[1]]), names(dfs[[2]])))

    # Generate group_clusterID dataset
    group_clusterID <- final_df %>%
      dplyr::select(contains("cluster")) %>%
      dplyr::select(all_of(
        names(.)[order(as.numeric(stringr::str_extract(names(.), "(?<=Group)\\d+")))]
      ))

    # Tidy final dataset by ordering and then removing cluster labels
    final_df <- final_df %>%
      dplyr::select(-contains("cluster")) %>%
      dplyr::select(
        Subject_ID, Time,  # fixed ID columns
        all_of(
          names(.) %>%
            setdiff(c("Subject_ID", "Time")) %>%
            .[order(as.numeric(str_extract(., "\\d+")))]  # sort numerically
        )
      )

      return(list(
        "Simulated Data" = final_df,
        "Cluster ID per participant per group" = group_clusterID,
        "Group ID per column" = group
      ))

    } else {
      # Single group: just return indiv_clust per participant
      sim_data <- as.data.frame(sim_data)
      return(list(
        "Simulated Data" = sim_data,
        "Cluster ID" = indiv_clust
      ))
    }
}