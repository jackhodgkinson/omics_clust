# simulateLCMM-test.R

# Load libraries
library(tidyverse)

# Load functions
source("simulateLCMM.R")

# Initialise parameters
n_groups <- 2
n_clust <- 2
n_col <- 10
n_indiv <- 419
timepoints <- c(-1.26, -0.345, 0.590, 1.53)
timepoints_sd  <- c(0.101, 0.0545, 0.0517, 0.0478)
seed <- 4881

# Cluster specific params 
params <- list(
  cluster1 = list(
    fixed_params = list(rnorm(n_col, 4.5, 0.2), rnorm(n_col, -0.1, 0.1), rnorm(n_col, -0.1, 0.1)),
    random_cov = matrix(c(mean(rnorm(n_col, 0.1, 0.05))^2, 0, 0, mean(rnorm(n_col, 0.02, 0.05))^2), 
                        nrow = 2, byrow = TRUE),
    resid_sd = rnorm(n_col, 0.5, 0.25)
  ),
  cluster2 = list(
    fixed_params = list(rnorm(n_col, -4.5, 0.4), rnorm(n_col, -0.1, 0.1), rnorm(n_col, -0.1, 0.1)),
    random_cov = matrix(c(mean(rnorm(n_col, 0.15, 0.6))^2, 0, 0, mean(rnorm(n_col, 0.02, 0.08))^2), 
                        nrow = 2, byrow = TRUE),
    resid_sd = rnorm(n_col, 0.3, 0.1)
  )
)

# ==== Simulating without Subjects ====
# Generate data
sim_data <- simulateLCMM(subject_data = NULL, ID = NULL, Time = NULL, timepoints, n_clust, n_groups, params, 
                        n_indiv, n_col, 4881, timepoint_noise = TRUE, timepoint_sd = timepoints_sd, cluster_labels = NULL, 
                        equal_clust = FALSE)

# Gather the data
sim_dat <- sim_data$`Simulated Data`
clusters <- sim_data$`Cluster ID per participant per group`
group <- sim_data$`Group ID`

# Plot the data
subjects <- unique(sim_dat$Subject_ID)
colnames(sim_dat)[-(1:2)] <- paste0("Protein", seq_len(n_col))

indiv_clust_df <- data.frame(Subject_ID = subjects, cluster = clusters$Group2_clusterID)

plot_data <- sim_dat %>%
  dplyr::select(Subject_ID, Time, everything()) %>%
  left_join(indiv_clust_df, by = "Subject_ID")

set.seed(seed)
sample <- sample(subjects, size = 50)
plot_data2 <- plot_data %>% filter(Subject_ID %in% sample)

ggplot(plot_data2, aes(x = Time, y = Protein1, group = Subject_ID)) +
  geom_line(color = "grey70", alpha = 0.4) +      # faint grey lines
  geom_point(aes(color = as.factor(cluster)), size = 2) +  # points coloured by cluster
  labs(title = "Longitudinal Plot by Cluster",
       x = "Time",
       y = "Protein Abundance",
       color = "Cluster Group") +
  theme_minimal()

# Produce heatmap - focused on group 2 and timepoint 1
sim_dat_filt <- sim_dat %>% filter(abs(Time - -1.26) <= 0.5)
sim_mat <- as.matrix(sim_dat_filt %>% dplyr::select(-Subject_ID, -Time))
annotationRow <- as.data.frame(clusters[["Group2_clusterID"]])
names(annotationRow) <- "Clusters"

annotationCol <- as.data.frame(group)
names(annotationCol) <- "ProteinClusters"
rownames(annotationCol) <- colnames(sim_mat)

rownames(sim_mat) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
annotationCol$ProteinClusters <- as.factor(annotationCol$ProteinClusters)
pheatmap::pheatmap(sim_mat[order(clusters[["Group2_clusterID"]]),
                                  order(annotationCol$ProteinClusters)], 
                   annotation_row = annotationRow,
                   annotation_col = annotationCol[order(annotationCol$ProteinClusters), , drop = FALSE], 
                   cluster_rows = F,
                   cluster_cols = F)

# Produce heatmap - focused on group 1 and timepoint 4
sim_dat_filt <- sim_dat %>% filter(abs(Time - 1.53) <= 0.5)
sim_mat <- as.matrix(sim_dat_filt %>% dplyr::select(-Subject_ID, -Time))
annotationRow <- as.data.frame(clusters[["Group1_clusterID"]])
names(annotationRow) <- "Clusters"

annotationCol <- as.data.frame(group)
names(annotationCol) <- "ProteinClusters"
rownames(annotationCol) <- colnames(sim_mat)

rownames(sim_mat) <- rownames(annotationRow)
annotationRow$Clusters <- as.factor(annotationRow$Clusters)
annotationCol$ProteinClusters <- as.factor(annotationCol$ProteinClusters)
pheatmap::pheatmap(sim_mat[order(clusters[["Group1_clusterID"]]),
                           order(annotationCol$ProteinClusters)], 
                   annotation_row = annotationRow,
                   annotation_col = annotationCol[order(annotationCol$ProteinClusters), , drop = FALSE], 
                   cluster_rows = F,
                   cluster_cols = F)

# ==== Simulating with Subjects ====
# Generate subjects
subj_data <- data.frame(SubjectID = rep(1:100, each = 4), 
                        Time = unlist(replicate(100, c(12, 20, 28, 36) + 
                                           rnorm(4, mean = 0, 
                                                 sd = c(0.842, 0.485, 0.419, 0.394)), 
                                         simplify = FALSE)),
                        GA = rep(c(12, 20, 28, 36), times = 4))

# Generate data
sim_data2 <- simulateLCMM(subject_data = subj_data, ID = NULL, Time = NULL, timepoints = NULL, n_clust, n_groups, params, 
                         100, n_col, 4881, cluster_labels = NULL, 
                         equal_clust = FALSE)

# Gather the data
sim_dat2 <- sim_data2$`Simulated Data`
clusters2 <- sim_data2$`Cluster ID per participant per group`
groups2 <- sim_data2$`Group ID`