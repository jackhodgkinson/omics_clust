## clusterofclusters.R
clusterofclusters <- function(moc,                                              # Matrix of Clusters - N X C data matrix, where C is the total number of clusters.
                              k = 1:9,                                          # Number of clusters. Default is to loop through k = 1 to k = 9.
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
                              random_seed = NULL                                # Set random seed for reproducibility. Default is NULL
                              ) 
  {   # Install relevant packages
      library(coca)
  
      # Set random seed 
      if(!is.null(random_seed)){
        set.seed(random_seed)
      }
      
      # Intialise output list
      output <- list()
      
      n <- dim(moc)[1]
      
      if (is.null(k)) {
        stop("Parameter 'k' is required. Please specify a value or range for 'k'.")
      }
      
      if (length(k) != 1 & choiceKmethod == "silhouette") {
        consensusMatrix <- array(NA, c(n, n, max(k) - 1))
        clLabels <- array(NA, c(max(k) - 1, n))
        
        for (i in seq_len(max(k) - 1) + 1) {
          ### Step 1. Compute the consensus matrix
          consensusMatrix[, , i - 1] <-
            consensusCluster(moc, i, B = N, pItem, clMethod = ccClMethod,
                             dist = ccDistHC, maxIterKM = max.iter)
          ### Step 2. Use hierarchical clustering on the consensus matrix
          distances <- stats::as.dist(1 - consensusMatrix[, , i - 1])
          hClustering <- stats::hclust(distances, method = hclustMethod)
          clLabels[i - 1, ] <- stats::cutree(hClustering, i)
        }
        
        K <- maximiseSilhouette(consensusMatrix, clLabels, max(k), savePNG,
                                fileName, widestGap = widestGap, dunns = dunns,
                                dunn2s = dunn2s)$K
        
      } else if (length(k) != 1 & choiceKmethod == "AUC") {
        consensusMatrix <- array(NA, c(n, n, max(k) - 1))
        areaUnderTheCurve <- rep(NA, max(k) - 1)
        
        for (i in seq_len(max(k) - 1) + 1) {
          ### Step 1. Compute the consensus matrix ###
          consensusMatrix[, , i - 1] <-
            consensusCluster(moc, i, B = N, pItem, clMethod = ccClMethod,
                             dist = ccDistHC, maxIterKM = max.iter)
          ### Step 2. Compute area under the curve ###
          areaUnderTheCurve[i - 1] <- computeAUC(consensusMatrix[, , i - 1])
        }
        
        K <- chooseKusingAUC(areaUnderTheCurve, savePNG, fileName)$K
        
      } else if (length(k) != 1) {
        stop("Method to choose number of clusters has not been recognised.
            Please make sure that it is either `silhouette` or `AUC`.")
      } else if (length(k) == 1) {
        consensusMatrix <- NULL
        K <- k
      }
      
      if (verbose)
        print(paste("K =", k, sep = " "))
      
      ### Step 1. Compute the consensus matrix ###
      if (!is.null(consensusMatrix)) {
        output$consensusMatrix <- consensusMatrix[, , K - 1]
      } else {
        output$consensusMatrix <- consensusCluster(moc, K, B = N, pItem,
                                                   clMethod = ccClMethod,
                                                   dist = ccDistHC,
                                                   maxIterKM = max.iter)
      }
      
      ### Step 2. Use hierarchical clustering on the consensus matrix ###
      distances <- stats::as.dist(1 - output$consensusMatrix)
      hClustering <- stats::hclust(distances, method = hclustMethod)
      output$clusterLabels <- stats::cutree(hClustering, K)
      
      output$K <- K
      
      if (returnAllMatrices)
        output$consensusMatrices <- consensusMatrix
      
      return(output)
      
    }