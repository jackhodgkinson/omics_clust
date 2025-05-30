# simulateLCMM
library(tidyverse)

n_groups <- 2
n_clust <- 3
N_col <- 25
n_indiv <- 448
timepoints <- c(12, 20, 28, 36)

# Cluster specific params 
params <- list(
  cluster1 = list(
    intercept = rnorm(N_col, 8, 1),
    sd_intercept = runif(N_col, 0.1, 0.5),
    sd_slope = runif(N_col, 0.02, 0.15),
    resid_sd = runif(N_col, 0.01, 0.1)
  ),
  cluster2 = list(
    intercept = rnorm(N_col, 0, 0.75),
    sd_intercept = runif(N_col, 0.15, 0.6),
    sd_slope = runif(N_col, 0.02, 0.08),
    resid_sd = runif(N_col, 0.05, 0.15)
  ),
  cluster3 = list(
    intercept = rnorm(N_col, -8, 0.6),
    sd_intercept = runif(N_col, 0.05, 0.4),
    sd_slope = runif(N_col, 0.04, 0.09),
    resid_sd = runif(N_col, 0.01, 0.17)
  )
)

# Test params
cluster_labels <- NULL 
equal_clust <- FALSE
seed <- 4881

# simulateLCMM <- function(subject_data = NULL,
#                          timepoints,
#                          n_clust,                                                 # Number of clusters
#                          n_groups,                                                # Number of groups of data, 1 by default
#                          cluster_params,                                          # List of distribution parameters per cluster
#                          n_indiv,                                                 # Number of individuals  
#                          n_col,                                                   # Number of columns in simulated data
#                          random_seed,                                             # Input random seed for reproducibility
#                          missing = FALSE,
#                          missing_perc = 0.05,
#                          missing_timepoints = NULL,                               # The index of timepoints you would like to add missing data to. Default is NULL, so selection is random.
#                          timepoint_perc = 0.1,                                    # The percentage of timepoints where you want missing observations removed 
#                          timepoint_noise = TRUE,                                  # Add noise to the timepoints. Default is TRUE.
#                          cluster_labels = NULL,                                   # Input cluster labels, NULL by default.
#                          group_labels = NULL,                                     # Input group labels, NA by default
#                          equal_clust = TRUE,                                      # If generating cluster labels, ensure each cluster has approx equal individuals
#                          equal_groups = TRUE,                                     # If n_groups > 1, ensure each group contains approx equal number of cols
#                          parallel_process = FALSE                                          # Use parallel processsing. Default is TRUE
# ){

  # Generate subject data if not provided
  #if (is.null(subject_data)) {
    
    ## Generate Subject / Time Data
    data_hlme <- data.frame(Subject_ID = rep(1:n_indiv, each = length(timepoints)),
                            Time = rep(timepoints, times = n_indiv))
    
    
    #if (missing == TRUE) {
      ## Remove some observations randomly
      # Get unique subject IDs
      set.seed(random_seed)
      unique_ids <- unique(data_hlme$Subject_ID)
      
      # Initialize empty index vector
      remove_idx <- c()
    
      # Randomly select subsets
      remove_1tp_ids <- sample(unique_ids, size = round(0.07 * n_indiv))
      remove_2tp_ids <- sample(setdiff(unique_ids, remove_1tp_ids), size = round(0.03 * n_indiv))
      
      
      
      # Remove 1 timepoint for some individuals
      for (id in remove_1tp_ids) {
        tp_to_remove <- sample(timepoints, 1)
        idx_to_remove <- which(data_hlme$Subject_ID == id & round(data_hlme$Time) == tp_to_remove)
        remove_idx <- c(remove_idx, idx_to_remove)
      }
      
      # Remove 2 timepoints for others
      for (id in remove_2tp_ids) {
        tps_to_remove <- sample(timepoints, 2)
        idx_to_remove <- which(data_hlme$Subject_ID == id & round(data_hlme$Time) %in% tps_to_remove)
        remove_idx <- c(remove_idx, idx_to_remove)
      }
      
      # Set selected rows to NA
      data_hlme$Time[remove_idx] <- NA
      
    #}
      
    # Identify non-missing indices
    non_missing_idx <- which(!is.na(data_hlme$Time))
    
    # Add noise
    data_hlme$Time[non_missing_idx] <- data_hlme$Time[non_missing_idx] +
      rnorm(length(non_missing_idx), mean = 0, sd = 0.25)
    
    # Generate cluster labels if not provided
    set.seed(random_seed)
    indiv_clust <- if (is.null(cluster_labels)){
      if (!equal_clust){
        sample(seq_len(n_clust), size = n_indiv, replace = TRUE, 
               prob = {p <- runif(n_clust); p / sum(p)})
      } else {
        sample(seq_len(n_clust), size = n_indiv, replace = TRUE)  
      } 
    } else cluster_labels

# Data simulation function
sim_longitud_data <- function (params, subject_data, ID, indiv_clust, timepoints) {
  
  sim_data <- matrix(NA, nrow = length(subject_data[[ID]]), ncol = N_col)
  colnames(sim_data) <- paste0("Protein", seq_len(N_col))
  
  for (k in seq_along(params)) {

    # Generate id for each cluster
    idx <- which(indiv_clust == k)
    if (length(idx) == 0) return(NULL)
    
    # Obtain parameters for cluster k
    param <- params[[paste0("cluster",k)]]
    
    for (i in idx) {
      
      # Get row indices in data_hlme for subject i
      subject_rows <- which(subject_data[[ID]] == i)
      if (length(subject_rows) == 0) next
    
      for (col_idx in seq_len(N_col)) {
        # Loop over each subject and time point
        baseline <- param$intercept[col_idx]+ rnorm(1, 0, param$sd_intercept[col_idx])
        indiv_slope <- rnorm(1, 0, param$sd_slope[col_idx])
        
        # Initial value
        current_value <- baseline
          
        for (t in seq_along(subject_rows)) {
          if (t > 1) {
            # Small step up/down from previous
            step_change <- rnorm(1, mean = indiv_slope, sd = param$resid_sd[col_idx])
            current_value <- current_value + step_change
          }
            
          # Add measurement noise
          sim_data[subject_rows[t], col_idx] <- rnorm(1, mean = current_value, sd = param$resid_sd[col_idx])
        }
      }
    }
  }
  sim_data <- cbind(subject_data, sim_data)
  rows_to_na <- which(rowSums(is.na(sim_data)) > 0)
  cols_to_na <- which(!(names(sim_data) %in% c(ID, "Time")))
  sim_data[rows_to_na, cols_to_na] <- NA
  return(sim_data)
}

sim_dat <- sim_longitud_data(params, data_hlme, "Subject_ID", indiv_clust, timepoints)


# Plot the data
subjects <- unique(sim_dat$Subject_ID)

indiv_clust_df <- data.frame(Subject_ID = subjects, cluster = indiv_clust)

plot_data <- sim_dat %>%
  left_join(indiv_clust_df, by = "Subject_ID") 

set.seed(seed)
sample <- sample(subjects, size = 50)
plot_data2 <- plot_data %>% filter(Subject_ID %in% sample)
  
ggplot(plot_data2, aes(x = Time, y = Protein3, group = Subject_ID)) +
  geom_line(color = "grey70", alpha = 0.4) +      # faint grey lines
  geom_point(aes(color = as.factor(cluster)), size = 2) +  # points colored by cluster
  labs(title = "Longitudinal Plot by Cluster",
       x = "Time",
       y = "Protein Abundance",
       color = "Cluster Group") +
  theme_minimal()
  