# 실루엣 계수 시각화 -------------------------------------------------------------------------
silhouette_visualize <- function(clustering_method = c("hclust","kmeans"), combine_similarity_fun, number_of_clusters, distance_type){
  if (clustering_method == "hclust"){
    hclust_result <- get_clustering_result(clustering_method, combine_similarity_fun, number_of_clusters, distance_type)
    best_alpha <- hclust_result[[2]]$alpha
    best_beta <- hclust_result[[2]]$beta
    best_gamma <- hclust_result[[2]]$gamma
    
    similarity <- combine_similarity(spatial_sim, time_sim, velocity_sim, par = c(best_alpha, best_beta, best_gamma))
    best_distance <- 1 - similarity
    sil <- silhouette(cutree(hclust_result[[3]], k = number_of_clusters), dist(as.matrix(best_distance)))
  }
  
  else if (clustering_method == "kmeans"){
    sil <- silhouette(best_kmeans$cluster, dist(as.matrix(best_dist)))
  }
  
  fviz_silhouette(sil)
}

## Usage Example -------------------------------------------------------------------------
combine_similarity_fun <- "weighted_average_similarity"
number_of_clusters <- 2
distance_type <- 1

silhouette_visualize("hclust", combine_similarity_fun, number_of_clusters, distance_type)