# optimal_hclust.R 
optimal_hclust <- function(dist_mat,
                           n_groups = NULL,
                           method, 
                           index) { 
  
  # Figure out optimal number of groups 
  set.seed(seed)
  
  # Try different linkage functions with different indices and see which is most accurate. Report this in thesis.
  results <- data.frame(Method = character(), 
                        Index = character(), 
                        NGroups = integer(), 
                        stringsAsFactors = FALSE)
  
  all_index_list <- list()
  
  for (m in method) {
    for (i in index) {
      try({
        opt <- NbClust(diss = dist_mat, distance = NULL,
                       method = m, min.nc = 2, max.nc = 9,
                       index = i)
        best_nc <- opt$Best.nc[[1]]
        results <- rbind(results, data.frame(Method = m,
                                             Index = i,
                                             BestNC = best_nc))
        index_df <- data.frame(NClust = as.numeric(names(opt$All.index)),
                               Value = as.numeric(opt$All.index),
                               Method = m,
                               Index = i)
        
        all_index_list[[paste(m, i, sep = "_")]] <- index_df
      }, silent = TRUE)
    }
  }
  
  # Join with true number of groups 
  opt_results <- cbind(results, n_groups)
  
  # Create data frame for plotting
  all_index_df <- bind_rows(all_index_list)
  
  # Plot graph 
  plot <- ggplot(all_index_df, aes(x = NClust, y = Value)) + 
    geom_line() + 
    facet_grid(Method ~ Index) + 
    labs(title = "NbClust Index Scores for Each Method and Index",
         x = "Number of Groups",
         y = "Index Value")
  
  print(plot)
  
  return(list(`Number of Clusters` = opt_results, 
              `Index Values` = all_index_df))
  }