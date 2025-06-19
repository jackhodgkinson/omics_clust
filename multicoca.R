# multicoca.R
multicoca <- function(moc_list,                                         # List of Matrix of Clusters
                      k = 2:9,                                          # Number of clusters. Default is to loop through k = 1 to k = 9.
                      N = 1000,                                         # Number of iterations of Consensus Clustering step.
                      max.iter = 1000,                                  # Maximum number of iterations for k-means clustering
                      pItem = 0.8,                                      # Proportion of items sampled at each iteration. 
                      hclustMethod = "average",                         # Agglomeration method to be used by the hclust function to perform hierarchical clustering on the consensus matrix. Can be "single","complete", "average", etc. For more details please see ?stats::hclust.
                      choiceKmethod = "silhouette",                     # Method used to choose the number of clusters if K is NULL, can be either "AUC" (area under the curve, work in progress) or "silhouette". Default is "silhouette".
                      ccClMethod = "kmeans",                            # Clustering method to be used by the Consensus Clustering algorithm (CC). Can be either "kmeans" for k-means clustering or "hclust" for hiearchical clustering. Default is "kmeans".
                      ccDistHC = "euclidean",                           # Distance to be used by the hiearchical clustering algorithm inside CC. Can be "pearson" (for 1 - Pearson correlation), "spearman" (for 1- Spearman correlation), or any of the distances provided in stats::dist() (i.e. "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"). Default is "euclidean".
                      savePNG = FALSE,                                  # Boolean. Save plots as PNG files. Default is FALSE.
                      fileName = "coca",                                # Boolean. If savePNG is TRUE, this is the string containing (the first part of) the name of the output files. Can be used to specify the folder path too. Default is "coca". The ".png" extension is automatically added to this string.
                      verbose = FALSE,                                  # Boolean. 
                      widestGap = FALSE,                                # Boolean. If TRUE, compute also widest gap index to choose best number of clusters. Default is FALSE.
                      dunns = FALSE,                                    # Boolean. If TRUE, compute also Dunn's index to choose best number of clusters. Default is FALSE.
                      dunn2s = FALSE,                                   # Boolean. If TRUE, compute also alternative Dunn's index to choose best number of clusters. Default is FALSE.
                      returnAllMatrices = FALSE,                        # Boolean. If TRUE, return consensus matrices for all considered values of K. Default is FALSE.
                      random_seed = NULL,                               # Set random seed for reproducibility. Default is NULL
                      parallel_process = FALSE                                   # Use parallel processing. Default is TRUE.
                      )  
{
  
  # Load necessary libraries
  library(coca)
  library(parallel)
  
  # Source function 
  source("clusterofclusters.R")
  source("numCores.R")
  
  # Set random seed 
  if(!is.null(random_seed)){
    set.seed(random_seed)
  }
  
  # Create an empty list to store the results
  results <- list()
  
  # Obtain number of available cores
  n_cores <- numCores() 
  
  # Run MultiCOCA with parallel processing
  if (parallel_process) {
    # Detect OS for parallel backend
    os_type <- .Platform$OS.type
    if (os_type != "windows") {
      results <- parallel::mclapply(1:length(moc_list), function(i) {
        result <- clusterofclusters(moc_list[[i]], 
                                    k = k,
                                    N = N, 
                                    max.iter = max.iter, 
                                    pItem = pItem, 
                                    hclustMethod = hclustMethod, 
                                    choiceKmethod = choiceKmethod, 
                                    ccClMethod = ccClMethod, 
                                    ccDistHC = ccDistHC, 
                                    savePNG = savePNG, 
                                    fileName = paste0(fileName, "_Group", i), 
                                    verbose = verbose, 
                                    widestGap = widestGap, 
                                    dunns = dunns, 
                                    dunn2s = dunn2s, 
                                    returnAllMatrices = returnAllMatrices,
                                    random_seed = random_seed,
                                    parallel_process = TRUE)
        
        # Name the results for each group
        return(result)
      })
    }
    else {
    
    # Create a cluster with available cores
    cl <- makeCluster(detectCores() - 1)
    
    # Export necessary variables and functions to the cluster
    clusterExport(cl, varlist = c("clusterofclusters"))
    
    # Run the parallelized operation using parLapply
    results <- parLapply(cl, 1:length(moc_list), function(i) {
        result <- clusterofclusters(moc_list[[i]], 
                                    k = k,
                                    N = N, 
                                    max.iter = max.iter, 
                                    pItem = pItem, 
                                    hclustMethod = hclustMethod, 
                                    choiceKmethod = choiceKmethod, 
                                    ccClMethod = ccClMethod, 
                                    ccDistHC = ccDistHC, 
                                    savePNG = savePNG, 
                                    fileName = paste0(fileName, "_Group", i), 
                                    verbose = verbose, 
                                    widestGap = widestGap, 
                                    dunns = dunns, 
                                    dunn2s = dunn2s, 
                                    returnAllMatrices = returnAllMatrices,
                                    random_seed = random_seed,
                                    parallel_process = TRUE)
        
        # Name the results for each group
        return(result)
    })
    
    # Stop the cluster after computation
    stopCluster(cl)
    
    # Name each group
    names(results) <- paste0("Group", 1:length(results))
    
    }
  }
  
  # Run MultiCOCA without parallel processing
  else {
    for (i in 1:length(moc_list)) {
      results[[i]] <- clusterofclusters(moc_list[[i]], 
                                        k = k, 
                                        N = N, 
                                        max.iter = max.iter, 
                                        pItem = pItem, 
                                        hclustMethod = hclustMethod, 
                                        choiceKmethod = choiceKmethod, 
                                        ccClMethod = ccClMethod, 
                                        ccDistHC = ccDistHC, 
                                        savePNG = savePNG, 
                                        fileName = paste0(fileName, "_Group", i), 
                                        verbose = verbose, 
                                        widestGap = widestGap, 
                                        dunns = dunns, 
                                        dunn2s = dunn2s, 
                                        returnAllMatrices = returnAllMatrices,
                                        random_seed = random_seed,
                                        parallel_process = FALSE)
      names(results)[[i]] <- paste0("Group", i)
    }
  }
  
  # Return the results list
  return(results)
}
