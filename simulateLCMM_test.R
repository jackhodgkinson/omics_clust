# simulateLCMM-test.R

# Load libraries
library(tidyverse)

# Load functions
source("simulateLCMM.R")

# Initialise parameters
n_groups <- 3
n_clust <- 3
N_col <- 25
n_indiv <- 448
timepoints <- c(12, 20, 28, 36)
seed <- 4881

# Cluster specific params 
params <- list(
  cluster1 = list(
    fixed_intercept = rnorm(N_col, 8, 1),
    fixed_slope = runif(N_col, -0.1, 0.1),
    random_intercept = runif(N_col, 0.1, 0.5),
    random_slope = runif(N_col, 0.02, 0.15),
    random_cov = matrix(c(mean(runif(N_col, 0.1, 0.5))^2, 0, 0, mean(runif(N_col, 0.02, 0.15))^2), 
                        nrow = 2, byrow = TRUE),
    resid_sd = runif(N_col, 0.1, 0.3)
  ),
  cluster2 = list(
    fixed_intercept = rnorm(N_col, 0, 0.75),
    fixed_slope = runif(N_col, -0.1, 0.1),
    random_intercept = runif(N_col, 0.15, 0.6),
    random_slope = runif(N_col, 0.02, 0.08),
    random_cov = matrix(c(mean(runif(N_col, 0.15, 0.6))^2, 0, 0, mean(runif(N_col, 0.02, 0.08))^2), 
                        nrow = 2, byrow = TRUE),
    resid_sd = runif(N_col, 0.15, 0.25)
  ),
  cluster3 = list(
    fixed_intercept = rnorm(N_col, -8, 0.6),
    fixed_slope = runif(N_col, -0.1, 0.1),
    random_intercept = runif(N_col, 0.05, 0.4),
    random_slope = runif(N_col, 0.04, 0.09),
    random_cov = matrix(c(mean(runif(N_col, 0.05, 0.4))^2, 0, 0, mean(runif(N_col, 0.04, 0.09))^2), 
                        nrow = 2, byrow = TRUE),
    resid_sd = runif(N_col, 0.21, 0.27)
  )
)

# Generate data
sim_data <- simulateLCMM(subject_data = NULL, timepoints, n_clust, n_groups, params, 
                        n_indiv, n_col = N_col, seed, missing = TRUE, timepoint_perc = 0.2, cluster_labels = NULL, 
                        equal_clust = FALSE)

# Plot the data
sim_dat <- sim_data$`Simulated Data`
clusters <- sim_data$`Cluster ID per participant per group`
group <- sim_data$`Group ID`

subjects <- unique(sim_dat$Subject_ID)
colnames(sim_dat)[-(1:2)] <- paste0("Protein", seq_len(N_col))

indiv_clust_df <- data.frame(Subject_ID = subjects, cluster = cluster)

plot_data <- sim_dat %>%
  left_join(indiv_clust_df, by = "Subject_ID") 

set.seed(seed)
sample <- sample(subjects, size = 50)
plot_data2 <- plot_data %>% filter(Subject_ID %in% sample)

ggplot(plot_data2, aes(x = Time, y = Protein16, group = Subject_ID)) +
  geom_line(color = "grey70", alpha = 0.4) +      # faint grey lines
  geom_point(aes(color = as.factor(cluster)), size = 2) +  # points coloured by cluster
  labs(title = "Longitudinal Plot by Cluster",
       x = "Time",
       y = "Protein Abundance",
       color = "Cluster Group") +
  theme_minimal()

# Produce heatmap - focused on group 1
sim_data_clust <- cbind(sim_dat, clusters[1])
sim_data_clust <- sim_data_clust %>%
  filter(abs(Time - 12) <= 1)
cluster <- sim_data_clust %>%
  select(Group1_clusterID) %>%
cluster <- as.numeric(cluster$Group1_clusterID)
sim_mat <- as.matrix(sim_data_clust %>% select(-Group1_clusterID, -Subject_ID, -Time))
annotationRow <- as.data.frame(cluster)
names(annotationRow) <- "Clusters"

annotationCol <- as.data.frame(as.factor(group))
names(annotationCol) <- "ProteinClusters"
rownames(annotationCol) <- colnames(sim_mat)

rownames(sim_mat) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
annotationCol$ProteinClusters <- as.factor(annotationCol$ProteinClusters)
pheatmap::pheatmap(sim_mat[order(cluster),
                                order(annotationCol$ProteinClusters)], 
                   annotation_row = annotationRow,
                   annotation_col = annotationCol[order(annotationCol$ProteinClusters), , drop = FALSE], 
                   cluster_rows = F)

