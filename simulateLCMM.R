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
                         # missing = FALSE,                                         # Boolean. Default is FALSE.
                         # missing_perc = 0.05,                                     # Percentage of missing data if missing == TRUE. Default is 0.05.
                         # missing_timepoints = NULL,                               # The index of timepoints you would like to add missing data to. Default is NULL, so selection is random.
                         # timepoint_perc = 0.1,                                    # The proportion of participants missing 1:n-1 timepoints. Can be a scalar or a vector of length n-1 timepoints. Default is 0.1
                         timepoint_noise = TRUE,                                  # Add noise to the timepoints. Default is TRUE. Only triggered if timepoint is not NULL.
                         timepoint_sd = 0,                                        # Specify the standard deviation of the timepoints to add noise. Can be a single number of a vector the same length as timepoints. Default is 0. Only triggered if timepoint is not NULL.
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
  # if (!is.logical(missing)) stop("`missing` must be TRUE or FALSE.")
  
  # if (!is.null(missing_timepoints)) {
  #   if (!all(missing_timepoints %in% timepoints)) {
  #     stop("All `missing_timepoints` must be in `timepoints`.")
  #   }
  # } else {
  #   if (!is.numeric(timepoint_perc) || any(timepoint_perc < 0)) {
  #     stop("timepoint_perc must be numeric and non-negative")
  #   }
  # }
  # 
  # if (length(missing_perc) != 1 && !is.null(missing_timepoints) &&
  #     length(missing_perc) != length(missing_timepoints)) {
  #   stop("`missing_perc` must be a single value or match length of `missing_timepoints`.")
  # }
  # 
  # if (any(missing_perc < 0 | missing_perc > 1)) {
  #   stop("`missing_perc` values must be between 0 and 1.")
  # }
  
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
  
  # ==== APPLY MISSINGNESS ====
  # if (missing) {
  #   set.seed(random_seed)
  #   subj <- unique(data_hlme$Subject_ID)
  #
  #   n_miss <- numeric(length(timepoint_perc))
  #
  #   # Calculate number of participants missing at least n-1,...,1 timepoint
  #   for (k in seq_along(timepoint_perc)) {
  #       n_miss[k] <- round(timepoint_perc[k] * length(subj))
  #   }
  #   n_miss <- c("0" = length(subj) - sum(n_miss), setNames(n_miss, seq_along(n_miss)))
  #
  #   # Initialise the target of missing counts per timepoint
  #   target_miss_tp <- round(missing_perc * length(subj))
  #   names(target_miss_tp) <- timepoints
  #
  #   # Extend n_miss to the number of participants and permute randomly.
  #   n_miss <- rep(as.integer(names(n_miss)), times = n_miss)
  #   set.seed(random_seed)
  #   n_miss <- sample(n_miss, length(n_miss), replace = FALSE)
  #
  #   # Matrix
  #   mat <- matrix(runif(length(n_miss * )))
  #
  #   # Matrix of missingness
  #   missing_mat <- matrix(FALSE, nrow = length(subj), ncol = length(timepoints))
  #
  #
  #   # Loop over each timepoint and mark n_miss timepoints as missing using target_miss_tp
  #   for (i in seq_along(subj)){
  #     if (n_miss[i] == 0) next
  #     tps <- names(which(colSums(missing_mat) < target_miss_tp))
  #
  #     n_assign <- min(n_miss[i], length(tps))
  #
  #     tp_assign <- sample(tps, n_assign)
  #
  #     missing_mat[i, tp_assign] <- TRUE
  #   }
  #
  #   return(missing_mat)
  # }
  #
  # if (missing) {
  #   set.seed(random_seed)
  #   all_subjects <- unique(data_hlme$Subject_ID)
  # 
  #   # Initialize rows to remove
  #   rows_to_remove <- integer(0)
  # 
  #   # Store missing selections as a data.frame: Subject_ID + Time
  #   missing_log <- data.frame(
  #     Subject_ID = character(0),
  #     Time = numeric(0)  # or character(0) if your timepoints are non-numeric
  #   )
  # 
  #   # 1. Select participants to have missing data at each timepoint according to missing_perc
  #   for (i in seq_along(timepoints)) {
  #     tp <- timepoints[i]
  #     perc <- missing_perc[i]
  # 
  #     tp_rows <- which(data_hlme$Time == tp)
  #     tp_subjects <- unique(data_hlme$Subject_ID[tp_rows])
  # 
  #     n_missing <- round(perc * length(tp_subjects))
  #     miss_subj <- sample(tp_subjects, size = n_missing, replace = FALSE)
  # 
  #     # Log each (Subject_ID, Time) pair
  #     missing_log <- rbind(missing_log, data.frame(
  #       Subject_ID = miss_subj,
  #       Time = tp
  #     ))
  # 
  #     tp_remove <- which(data_hlme$Subject_ID %in% miss_subj & data_hlme$Time == tp)
  #     rows_to_remove <- c(rows_to_remove, tp_remove)
  #   }
  # 
  #   # Correct number of participants with multiple missing timepoints to match timepoint_perc
  #   subject_missing_counts <- table(missing_log$Subject_ID)
  #   multi_missing <- names(subject_missing_counts[subject_missing_counts > 1])
  #   n_total_multi <- length(multi_missing)
  #   n_target_multi <- round(timepoint_perc * length(subject_missing_counts))
  # 
  #   # Randomly select a subset to keep
  #   if (n_total_multi > n_target_multi) {
  #     keep <- sample(multi_missing, size = n_target_multi)
  #     drop <- setdiff(multi_missing, keep)
  #   }
  # 
  #   keep_rows <- data.frame(Subject_ID = character(), Time = numeric(), stringsAsFactors = FALSE)
  # 
  #   for (id in drop){
  #     subj_log <- subset(missing_log, Subject_ID == id)
  #     keep_tp <- subj_log[sample(1:nrow(subj_log), size = 1),]
  #     keep_rows <- rbind(keep_rows, keep_tp)
  #   }
  # 
  #   # Correct missing_log
  #   missing_log <- rbind(missing_log[!(missing_log$Subject_ID %in% drop), ], keep_rows)
  # 
  #   # Correct number of participants with missing data at each timepoint by selecting new participants
  #   for (i in seq_along(timepoints)) {
  #     tp <- timepoints[i]
  #     target_miss <- round(missing_perc[i] * length(unique(data_hlme$Subject_ID[data_hlme$Time == tp])))
  # 
  #     # Current missing subjects at timepoint
  #     curr_miss <- unique(missing_log$Subject_ID[missing_log$Time == tp])
  #     n_curr_miss <- length(curr_miss)
  # 
  #     # Number of new subjects to add
  #     n_add <- target_miss - n_curr_miss
  # 
  #     # Eligble participants for new missing at this tp
  #     if (n_add > 0) {
  #       curr_miss2 <- unique(missing_log$Subject_ID)
  #       av_subj <- setdiff(unique(data_hlme$Subject_ID[data_hlme$Time == tp]), curr_miss2)
  # 
  #       if (length(av_subj) >= n_add) {
  #         new_miss <- sample(av_subj, n_add, replace = FALSE)
  #       } else {
  #         new_miss <- av_subj
  #       }
  # 
  #       missing_log <- rbind(missing_log, data.frame(Subject_ID = new_miss, Time = tp))
  # 
  #       new_rows <- which(data_hlme$Subject_ID %in% new_miss & data_hlme$Time == tp)
  #       rows_to_remove <- c(rows_to_remove, new_rows)
  #     }
  #   }
  # 
  #   ## TEST
  # 
  #   # # Get full list of timepoints
  #   all_timepoints <- sort(unique(data_hlme$Time))
  # 
  #   # Total number of participants (assuming all timepoints originally present)
  #   n_participants <- length(all_subjects)
  # 
  #   # 1. % of participants with missing data at each timepoint
  #   missing_by_tp <- sapply(all_timepoints, function(tp) {
  #     missing_subjects <- unique(missing_log$Subject_ID[missing_log$Time == tp])
  #     round(100 * length(missing_subjects) / n_participants, 2)
  #   })
  #   names(missing_by_tp) <- all_timepoints
  # 
  #   # # 2. % of participants with multiple timepoints missing
  #   subject_missing_counts <- table(missing_log$Subject_ID)
  #   n_missing_any <- length(subject_missing_counts)
  #   n_missing_multi <- sum(subject_missing_counts > 1)
  #   multi_tp_missing_perc <- round(100 * n_missing_multi / n_missing_any, 2)
  #   return(multi_tp_missing_perc)
  
  
  # # 3. Remove multiple timepoints from selected individuals
  # for (subj in multi_missing_participants) {
  #   subj_rows <- which(data_hlme$Subject_ID == subj)
  #   subj_tp <- unique(data_hlme$Time[subj_rows])
  #
  #   if (length(subj_tp) > 2) {
  #     n_remove <- sample(2:(n_tp - 1), 1)
  #   } else {
  #     n_remove <- 2
  #   }
  #
  #   tp_to_remove <- tail(subj_rows, n_remove)
  #   rows_to_remove <- c(rows_to_remove, tp_to_remove)
  # }
  #
  #
  # for (subj in single_missing_participants) {
  #   subj_rows <- which(data_hlme$Subject_ID == subj)
  #   tp_to_remove <- tail(subj_rows, 1)
  #   rows_to_remove <- c(rows_to_remove, tp_to_remove)
  # }
  #
  #
  #
  # # Drop the rows from data_hlme to remove missing timepoints
  # data_hlme <- data_hlme[-rows_to_remove, ]
  # }
  
  
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
  #
  
  # ==== SIMULATION FUNCTION ====
  
  n_cores <- numCores()
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
          
          re_contrib <- b_i[1]  # intercept
          
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
      view <- cbind(sim_data[, c(ID, TimeVar), drop = FALSE], view)
      name <- paste0("Group", g)
      group_list[[name]] <- view
    }
    
    # Empty list for the permuted groups
    group_permute <- list()
    
    # Permute the columns per group
    for (i in seq_along(names(group_list))) {
      name <- names(group_list)[i]
      group_df <- group_list[[name]]
      
      permute_col <- function(data, ID, TimeVar, cluster) {
        clust_data <- as.data.frame(cbind(unique(data[[ID]]), cluster))
        colnames(clust_data)[1] <- ID
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
    
      group_permute[[name]] <- permute_col(group_df, ID, TimeVar, indiv_clust)
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
      df[order(df[[ID]], df[[TimeVar]]), ]
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