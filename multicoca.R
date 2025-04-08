# multicoca.R
multicoca <- function(data,              # Input as data frame or a list of data frames
                      random_seed,       # Input random seed for reproducibility
                      N = 2000,          # Number of iterations of Consensus Clustering step
                      max.iter = 2000    # Maximum number of iterations for k-means clustering
                      )
{
  
  # Load packages
  library(lcmm)
  library(mclust)
  library(NbClust)
  
  # Set seed 
  set.seed(random_seed)
  
  # Fit LCMM model
  lcmm_mod <- hlme(~ - 1.)
  
  # Create similarity matrix
  sim_mat <- matrix(NA, ncol = ncol(data), nrow = nrow(data))
  for (i in 1:ncol(data)){
    for(j in 1:ncol(data)){
      sim_mat[i,j] <- adjustedRandIndex()   
    }
  }
  
  # Convert to dissimiliarity matrix 
  dissim <- 1 - sim_mat
  
  # Obtain optimal number of clusters
  opt_clust <- NbClust()
  
  # Create dendrogram and split it 
  hier_clust <- hclust()
  
  ## This is where we cluster the proteins and then feed each group into COCA. 
  
  # Run COCA 
  
  # Output: list of n_groups dataframes with cluster labels. 
  
  
}