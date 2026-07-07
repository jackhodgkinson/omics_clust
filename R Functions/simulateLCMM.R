# simulateLCMM
# NOTE: If participants want missing data, they need to provide this through subject_data
simulateLCMM <- function(subject_data = NULL,                                     # Provide subject data. Default is NULL and function will generate subject data with no missingness.
                         ID = NULL,                                               # If specifying subject_data, provide ID column. Default is NULL.
                         Time = NULL,                                             # If specifying subject_data, provide Time column. Default is NULL
                         timepoints = NULL,                                       # Timepoints for longitudinal data if not provided in subject_data. Default is NULL.
                         n_clust,                                                 # Number of clusters
                         n_groups,                                                # Number of groups of data, 1 by default
                         cluster_params,                                          # List of distribution parameters per cluster
                         n_indiv,                                                 # Number of individuals. 
                         n_col,                                                   # Number of columns in simulated data
                         random_seed,                                             # Input random seed for reproducibility
                         timepoint_noise = TRUE,                                  # Add noise to the timepoints. Default is TRUE. Only triggered if timepoint is not NULL.
                         timepoint_sd = 0,                                        # Specify the standard deviation of the timepoints to add noise. Can be a single number of a vector the same length as timepoints. Default is 0. Only triggered if timepoint is not NULL.
                         cluster_labels = NULL,                                   # Input cluster labels, NULL by default.
                         group_labels = NULL,                                     # Input group labels, NULL by default
                         equal_clust = TRUE,                                      # Boolean. If generating cluster labels, ensure each cluster has approx equal individuals. Default is TRUE.
                         equal_groups = TRUE,                                     # Boolean. If n_groups > 1, ensure each group contains approx equal number of columns. Default is TRUE.
                         parallel_process = TRUE                                  # Boolean. Use parallel processsing. Default is TRUE
){
  
  # ==== SET SEED FOR REPRODUCIBILITY ==== 
  RNGkind("L'Ecuyer-CMRG")
  set.seed(random_seed)
  
  # ==== PARAMETER VALIDATION ====
  if(length(cluster_params) != n_clust) {
    stop("Length of cluster parameters must match n_clust")
  }
  if (is.null(timepoints)) {
    if (is.null(subject_data)) {
      stop("`timepoints` must be specified if `subject_data` is NULL.")
    }
  } else if (!is.numeric(timepoints) || length(timepoints) < 1) {
    stop("`timepoints` must be a non-empty numeric vector.")
  }  
  
  if (!is.null(subject_data)) {
    
    # Handle ID column
    if (is.null(ID)) {
      id_col <- grep("ID", colnames(subject_data), value = TRUE, ignore.case = TRUE)
      
      if (length(id_col) != 1) {
        stop("`subject_data` must contain exactly one column with 'ID' in its name or you must specify the `ID` argument explicitly.")
      }
      
      ID <- id_col
      
    } else {
      if (!ID %in% colnames(subject_data)) {
        stop(paste0("`ID` column '", ID, "' not found in subject_data."))
      }
    }
    
    # Handle Time column
    if (is.null(Time)) {
      # No Time specified, try to find exactly one Time column
      time_cols <- grep("Time", colnames(subject_data), value = TRUE, ignore.case = TRUE)
      
      if (length(time_cols) == 0) {
        stop("`subject_data` must contain at least one column with 'Time' in its name or you must specify the `Time` argument explicitly.")
        
      } else if (length(time_cols) == 1) {
        TimeVar <- time_cols
        
      } else {
        stop(
          paste0(
            "Multiple columns with 'Time' found in `subject_data`: ",
            paste(time_cols, collapse = ", "),
            ". Please specify the exact `Time` column name using the `Time` argument."
          )
        )
      }
      
    } else {
      # User specified a Time column, check it exists
      if (!Time %in% colnames(subject_data)) {
        stop(paste0("`Time` column '", Time, "' not found in subject_data."))
      }
      TimeVar <- Time 
    }
    
    # Validate number of unique IDs matches n_indiv
    n_unique_ids <- length(unique(subject_data[[ID]]))
    
    if (n_unique_ids != n_indiv) {
      stop(
        paste0(
          "`n_indiv` (", n_indiv,
          ") does not match number of unique IDs in `subject_data` (", n_unique_ids, ")."
        )
      )
    }
    
  } else {
    # subject_data is NULL, handle accordingly
    if (is.null(timepoints)) {
      stop("`timepoints` must be specified if `subject_data` is NULL.")
    }
  }
  
  
  if (!is.numeric(n_clust) || n_clust < 1) stop("`n_clust` must be a positive integer.")
  if (!is.numeric(n_groups) || n_groups < 1) stop("`n_groups` must be a positive integer.")
  if (!is.list(cluster_params)) stop("`cluster_params` must be a list.")
  if (!is.numeric(n_indiv) || n_indiv < 1) stop("`n_indiv` must be a positive integer.")
  if (!is.numeric(n_col) || n_col < 1) stop("`n_col` must be a positive integer.")
  if (!is.numeric(random_seed)) stop("`random_seed` must be numeric.")

  # ==== GENERATE SUBJECT DATA ====  
  if (is.null(subject_data)) {
    # No subject_data provided: create default data frame with user-supplied timepoints and IDs
    data_hlme <- data.frame(
      Subject_ID = rep(1:n_indiv, each = length(timepoints)),
      Time = rep(timepoints, times = n_indiv)
    )
    ID <- "Subject_ID"
    TimeVar <- "Time"
  } else {
    
    
    # Check ID column presence
    if (!ID %in% colnames(subject_data)) {
      stop(paste0("`ID` column '", ID, "' not found in subject_data."))
    }
    
    # Check Time column presence
    if (!TimeVar %in% colnames(subject_data)) {
      stop(paste0("`Time` column '", TimeVar, "' not found in subject_data."))
    }
    
    data_hlme <- subject_data
  }
  
  # ==== ADD TIMEPOINT NOISE ====
  non_missing_idx <- which(!is.na(data_hlme$Time))
  
  if (!is.null(timepoints)) {
    if (timepoint_noise) {
      non_missing_idx <- which(!is.na(data_hlme$Time))
      
      if (length(timepoint_sd) == 1) {
        sd_vals <- rep(timepoint_sd, length(non_missing_idx))
      } else if (length(timepoint_sd) == length(timepoints)) {
        sd_vals <- timepoint_sd
      }
      
      data_hlme$Time[non_missing_idx] <- data_hlme$Time[non_missing_idx] +
        rnorm(length(non_missing_idx), mean = 0, sd = sd_vals)
    }
  }
  
  # ==== CLUSTER LABELS ====
  indiv_clust <- if (is.null(cluster_labels)) {
    if (!equal_clust) {
      sample(seq_len(n_clust), size = n_indiv, replace = TRUE,
             prob = {p <- runif(n_clust); p / sum(p)})
    } else {
      rep_len(rep(1:n_clust, length.out = n_indiv), n_indiv)
    }
  } else cluster_labels
  
  indiv_clust_long <- rep(NA_integer_, nrow(data_hlme))
  subject_ids <- data_hlme[[ID]]
  for (i in seq_along(indiv_clust)) {
    subj_id <- i  
    rows_subj <- which(subject_ids == subj_id)
    indiv_clust_long[rows_subj] <- indiv_clust[i]
  }
  
  # ==== SIMULATION FUNCTION ====
  
  n_cores <- parallel::detectCores()-2
  sim_longitud_data <- function(params, subject_data, ID, TimeVar, indiv_clust, n_col, parallel_proc = parallel_process) {
    
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
      
      # Validate fixed_params (a list of vectors)
      if (!("fixed_params" %in% names(p)) || !is.list(p$fixed_params)) {
        stop(paste(cluster_name, ": fixed_params must be a list"))
      }
      
      n_params <- length(p$fixed_params)
      if (n_params == 0) {
        stop(paste(cluster_name, ": fixed_params cannot be empty"))
      }
      
      # Check each fixed_params element length = n_col
      for (i in seq_len(n_params)) {
        if (length(p$fixed_params[[i]]) != n_col) {
          stop(paste(cluster_name, ": fixed_params element", i, "must be vector of length", n_col))
        }
      }
      
      # Validate resid_sd vector
      if (!("resid_sd" %in% names(p)) || length(p$resid_sd) != n_col) {
        stop(paste(cluster_name, ": resid_sd must be vector of length", n_col))
      }
      if (any(p$resid_sd < 0)) {
        stop(paste(cluster_name, ": resid_sd values must be non-negative"))
      }
      
      # Validate random_cov matrix
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
      if (any(eigenvalues < 0)) {
        stop(paste(cluster_name, ": random_cov must be positive semi-definite"))
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
    
    sim_cluster <- function(k, subject_data_local, params_local, ID_local, TimeVar_local, indiv_clust_local, n_col_local) {
      
      # Get parameters for the cluster and see how many fixed and random effects components 
      param <- params_local[[paste0("cluster", k)]]
      n_re <- nrow(param$random_cov)
      n_fe <- length(param$fixed_params)
      
      # See which subjects are being simulated based on cluster ID 
      cluster_subjects <- names(indiv_clust_local)[indiv_clust_local == k]
      if (length(cluster_subjects) == 0) return(NULL)
      
      # Initialise simulation matrix 
      sim_mat <- matrix(NA, nrow = nrow(subject_data_local), ncol = n_col_local)
      
      # Loop over each subject in the cluster
      for (subj_id in cluster_subjects) {
        
        # Gather information for each of the subjects in teh cluster
        subject_rows <- which((subject_data_local[[ID_local]] == subj_id) & (!is.na(subject_data_local[[TimeVar_local]])))
        if (length(subject_rows) == 0) next
        
        # Specify time variables
        time_vals <- subject_data_local[[TimeVar_local]][subject_rows]
        
        # Drawn random effects vector for this subject
        b_i <- MASS::mvrnorm(1, mu = rep(0, n_re), Sigma = param$random_cov)
        
        # Simulate outcome for each column
        outcome <- sapply(seq_len(n_col_local), function(col_idx) {
          
          re_contrib <- b_i[1]  
          
          if (n_re >= 2) {
            re_contrib <- re_contrib + b_i[2] * time_vals  # slope
          }
          if (n_re >= 3) {
            re_contrib <- re_contrib + b_i[3] * time_vals^2  # quadratic term
          }
          
          if (n_fe >= 1) {
            fe_contrib <- param$fixed_params[[1]][col_idx]
          }
          if (n_fe >= 2) {
            fe_contrib <- param$fixed_params[[1]][col_idx] + param$fixed_params[[2]][col_idx]*time_vals
          }
          if (n_fe >= 3) {
            fe_contrib <- param$fixed_params[[1]][col_idx] + param$fixed_params[[2]][col_idx]*time_vals + param$fixed_params[[3]][col_idx]*time_vals^2
          }
          
          resid_noise <- rnorm(1, 0, param$resid_sd)
          
          fe_contrib + re_contrib + resid_noise
          
        })
        
        sim_mat[subject_rows, ] <- outcome
      }
      
      sim_mat
    }
    
    if (parallel_proc) {
      # Detect OS for parallel backend
      os_type <- .Platform$OS.type
      clus <- NULL
      
      if (os_type == "windows") {
        clus <- parallel::makeCluster(n_cores)
        parallel::clusterSetRNGStream(clus, iseed = random_seed)  # sets streams per worker
        parallel::clusterExport(clus, varlist = c("params", "subject_data", "ID", "TimeVar", "indiv_clust", "n_col", "sim_cluster"), envir = environment())
        parallel::clusterEvalQ(clus, library(MASS))
        sim_results <- parallel::parLapply(
          clus,
          seq_along(params),
          function(k) sim_cluster(k, subject_data, params, ID, TimeVar, indiv_clust, n_col)
        )
        parallel::stopCluster(clus)
      } else {
        # Unix/macOS: use mclapply
        sim_results <- parallel::mclapply(
          seq_along(params),
          function(k) sim_cluster(k, subject_data, params, ID, TimeVar, indiv_clust, n_col),
          mc.cores = n_cores
        )
      }
    } else {
      sim_results <- lapply(seq_along(params), function(k) sim_cluster(k, subject_data, params, ID, TimeVar, indiv_clust, n_col))
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
  sim_data <- sim_longitud_data(cluster_params, data_hlme, ID, TimeVar,
                                indiv_clust, n_col, parallel_proc = parallel_process)
  
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
        group <- sample(seq_len(n_groups), n_col, replace = TRUE, prob = p)
        if (length(unique(group)) == n_groups) break
      }
    } else {
      group <- sample(seq_len(n_groups), n_col, replace = TRUE)
    }
    
    # Split the dataset into the different groups
    group_list <- list()
    split_data <- sim_data %>%
      dplyr::select(where(is.numeric), -any_of(names(data_hlme)))
    
    unique_groups <- unique(group)
    for (g in unique_groups[-length(unique_groups)]) {  # Exclude last group
      view <- split_data[, group == g, drop = FALSE]
      view <- cbind(sim_data[, c(ID, Time)], view)
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
        clust_data <- as.data.frame(cbind(cluster, unique(data[[ID]])))
        colnames(clust_data)[2] <- ID
        data_clust <- data %>%
          left_join(clust_data, by = ID)
        
        subject_counts <- data_clust %>%
          group_by(.data[[ID]]) %>%
          summarise(n_timepoints = n())
        
        data_clust_counts <- data_clust %>%
          left_join(subject_counts, by = ID)
        
        data_by_timepoint <- split(data_clust_counts, data_clust_counts$n_timepoints)
        names(data_by_timepoint) <- paste(seq_along(data_by_timepoint),"Timepoint(s)")
        data_by_timepoint <- lapply(data_by_timepoint, function(df) {df$n_timepoints <- NULL
        return(df)
        })
        
        data_by_timepoint <- lapply(data_by_timepoint, function(df) {
          
          orig_ID <- df[[ID]]
          orig_time <- df[[TimeVar]]
          permute_ID <- sample(unique(df[[ID]]))
          df <- df %>%
            arrange(match(df[[ID]], permute_ID), TimeVar) %>%
            mutate(!!ID := orig_ID,
                   !!TimeVar := orig_time)
          df
        })
        
        data2 <- bind_rows(data_by_timepoint) %>%
          arrange(ID, TimeVar)
        
        return(data2)
      }
      
      group_permute[[name]] <- permute_col(group_df, ID, indiv_clust, seed)
      names(group_permute[[name]])[names(group_permute[[name]]) == "cluster"] <- paste0(name, "_clusterID")
    }
    
    # Obtain data for non-permuted group
    non_permute_group_num <- unique_groups[length(unique_groups)]
    non_permute_view <- split_data[, group == non_permute_group_num, drop = FALSE]
    non_permute_group <- cbind(sim_data[, c(ID, TimeVar)], non_permute_view)
    non_perm_clust_data <- data.frame(
      cluster = indiv_clust
    )
    non_perm_clust_data[[ID]] <- unique(sim_data[[ID]])
    non_perm_data_clust <- non_permute_group %>%
      left_join(non_perm_clust_data, by = ID)
    names(non_perm_data_clust)[names(non_perm_data_clust) == "cluster"] <- 
      paste0("Group",non_permute_group_num, "_clusterID")
    
    # List of dataframes
    name <- paste0("Group", non_permute_group_num)
    named_non_perm <- setNames(list(non_perm_data_clust), name)
    dfs <- c(named_non_perm, group_permute)
    dfs <- lapply(dfs, function(df) {
      df[order(df[[ID]], df$Time), ]
    })
    
    # Join all dataframes avoiding duplication of ID and Time
    final_df <- reduce(dfs, full_join, by = intersect(names(dfs[[1]]), names(dfs[[2]])))
    
    # Generate group_clusterID dataset
    group_clusterID <- final_df %>%
      dplyr::select(all_of(ID), contains("cluster")) %>%
      dplyr::select(all_of(
        names(.)[order(as.numeric(stringr::str_extract(names(.), "(?<=Group)\\d+")))]
      )) %>%
      dplyr::distinct() %>%
      dplyr::select(-all_of(ID))
    
    # Tidy final dataset by ordering and then removing cluster labels
    final_df <- final_df %>%
      dplyr::select(-contains("cluster")) %>%
      dplyr::select(
        all_of(ID), TimeVar,  # fixed ID columns
        all_of(
          names(.) %>%
            setdiff(c(ID, "Time")) %>%
            .[order(as.numeric(str_extract(., "\\d+")))]  # sort numerically
        )
      )
    
    # Add back any columns from subject_data that were lost
    if (!is.null(subject_data)) {
      extra_cols <- setdiff(colnames(subject_data), colnames(final_df))
      if (length(extra_cols) > 0) {
        final_df <- dplyr::left_join(final_df, subject_data[, c(ID, TimeVar, extra_cols)], by = c(ID, TimeVar))
        final_df <- final_df %>%
          dplyr::select(
            which(!grepl("\\d", colnames(final_df))),
            which(grepl("\\d", colnames(final_df)))
          )
        
      }
    }
    
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