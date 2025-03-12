setwd("path/to/your/wd")

load_packages <- function(){
  #install.packages("fastcluster")  
  library(fastcluster) #hclust
  library(cluster) #실루엣 계수
  library(factoextra) #실루엣 계수 시각화
}

packageVersion("fastcluster")
packageVersion("cluster")
packageVersion("factoextra")


load_packages()

dropbox_data_folder_path <- "path/to/dropbox/folder"
#similarity 저장되어있는 폴더

data_folder_path <- "path/to/your/wd"
#결과 저장 폴더위치

sym_spatial <- read.csv(paste0(dropbox_data_folder_path, "result/2024.10.3/sym_spatial.csv"))
sym_time <- read.csv(paste0(dropbox_data_folder_path, "result/2024.10.3/sym_time.csv"))
velocity_sim <- read.csv(paste0(dropbox_data_folder_path, "result/2024.09.21/velocity_sim_939.csv"))
spatial_sim <- subset(sym_spatial, select = -X)
time_sim <- subset(sym_time, select = -X)
velocity_sim <- subset(velocity_sim, select = -X)
rm(sym_time)
rm(sym_spatial)


# combine_similarities ----------------------------------------------------
weighted_average_similarity <- function(spatial_sim=NULL, time_sim=NULL, velocity_sim=NULL,
                                        par=c(1,1,1)){
  alpha <- par[1]; beta <- par[2]; gamm <- par[3]
  return((as.matrix(spatial_sim)*alpha + as.matrix(time_sim)*beta + as.matrix(velocity_sim)*gamm)/(alpha + beta + gamm))
}



# alpha, beta, gamma
seq_vals <- seq(0, 1, by = 0.05)
combinations <- expand.grid(alpha = seq_vals, beta = seq_vals, gamma = seq_vals)
valid_combinations <- combinations[rowSums(combinations) == 1, ]
print(valid_combinations)


# grid seacrh -------------------------------------------------------------
clustering_result_folder_path <- paste0(data_folder_path,"클러스터링 결과 8/")

clustering <- function(cluster_model, combine_similarity_fun, number_of_clusters) {
  for (k in 2:number_of_clusters){
    best_silhouette <- 0
    best_params <- list()
    
    print(k)
    for (i in seq(1,nrow(valid_combinations))) {
      alpha <- valid_combinations[i,1]
      beta <- valid_combinations[i,2]
      gamma <- valid_combinations[i,3]
      
      print(paste(i, alpha, beta, gamma))
      # 중간 확인용 프린트
      
      similarity <- get(combine_similarity_fun)(spatial_sim, time_sim, velocity_sim, par = c(alpha, beta, gamma))
      
      ## Set distance ------------------------------------------------------------
      distance_matrix <- 1 - similarity; dist_type <- 1
      
      ## Kmeans ------------------------------------------------------------------
      # set.seed(123)
      # clustering_result <- kmeans(distance_matrix, centers = 3, nstart = 25)
      # sil <- silhouette(clustering_result$cluster, dist(distance_matrix))
      # avg_silhouette <- mean(sil[, 3])
      
      ## Hclust ------------------------------------------------------------------
      if (cluster_model == "hclust"){
        distance_matrix <- as.dist(distance_matrix)  
        clustering_result <- hclust(distance_matrix, method = "ward.D2")
        # k <- 3 #군집의수 
        clusters <- cutree(clustering_result, k = k)
        sil <- silhouette(clusters, distance_matrix)
        avg_silhouette <- mean(sil[, 3])
      }
      
      # 최적의 실루엣 계수 및 파라미터 저장
      if (avg_silhouette > best_silhouette) {
        best_silhouette <- avg_silhouette
        best_params <- list(alpha = alpha, beta = beta, gamma = gamma)
        best_clustering <- clustering_result
        best_dist <- distance_matrix
      }
      print(avg_silhouette)
    }
    
    hclust_clusters_results <- list(
      "best_silhouette" <- best_silhouette,
      "best_parameters" <- best_params,
      "best_clustering" <- best_clustering
    )
    
    rds_name <- paste(cluster_model, combine_similarity_fun, k, dist_type, "results.rds",sep="_")
    file_path <- paste(clustering_result_folder_path, rds_name, sep = "")
    saveRDS(hclust_clusters_results, file_path)
    
  }
}


## Run Grid Search ---------------------------------------------------------
seq_vals <- seq(0, 1, by = 0.1)
combinations <- expand.grid(alpha = seq_vals, beta = seq_vals, gamma = seq_vals)
valid_combinations <- combinations[rowSums(combinations) == 1, ]
print(valid_combinations)

set.seed(42)

clustering("hclust","weighted_average_similarity",10)

file_names <- list.files(clustering_result_folder_path)
for (file in file_names){
  print(file)
  result <- readRDS(paste0(clustering_result_folder_path,"/",file))
  print(result[[2]])
}

