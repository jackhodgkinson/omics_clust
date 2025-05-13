# groupARI.R
groupARI <- function(output,
                     true_clusters){
  group_names <- names(results)
  ari_values <- data.frame(Group = character(),
                           ARI = numeric(),
                           stringsAsFactors = FALSE
                           )
  for (group in group_names){
    ari_values <- rbind(ari_values, data.frame(Group = group,
                                               ARI = adjustedRandIndex(
                                                 true_clusters[[paste0(tolower(group),"_clusterid")]],
                                                 output[[group]]$clusterLabels)
    ))
  }
  
  return(ari_values)
}