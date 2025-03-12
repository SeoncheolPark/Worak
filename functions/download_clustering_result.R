library(fastcluster) #hclust
library(cluster) #실루엣 계수
library(factoextra) #실루엣 계수 시각화

# 드랍박스 데이터 저장 위치
dropbox_data_folder_path <- "set/dropbox/data/folder/"

# 로컬 데이터 저장 위치
data_folder_path <- "set/local/data/folder/"

# 로컬 클러스터링 결과 저장 위치
cluster_results_folder_path <- "set/local/clustering/result/folder/"

# 계산된 유사도 불러오기
spatial_sim <- read.csv(paste0(dropbox_data_folder_path, "result/2024.10.3/sym_spatial.csv"))
time_sim <- read.csv(paste0(dropbox_data_folder_path, "result/2024.10.3/sym_time.csv"))
velocity_sim <- read.csv(paste0(dropbox_data_folder_path, "result/2024.09.21/velocity_sim_939.csv"))
spatial_sim <- subset(spatial_sim, select = -X)
time_sim <- subset(time_sim, select = -X)
velocity_sim <- subset(velocity_sim, select = -X)

# 클러스터링 결과 불러오는 함수 -------------------------------------------------------------------------

get_clustering_result <- function(clustering_method = c("hclust","kmeans"), combine_similarity_fun, number_of_clusters, distance_type){
  rds_name <- paste(clustering_method, combine_similarity_fun, number_of_clusters, distance_type, "results.rds", sep = "_")
  file_path <- paste(cluster_results_folder_path,rds_name, sep = "")
  cluster_result <- readRDS(file_path)
  
  return(cluster_result)
}

## Usage Example -------------------------------------------------------------------------
combine_similarity_fun <- "weighted_average_similarity"
number_of_clusters <- 2
distance_type <- 1

# get clustering result
hclust_result <- get_clustering_result("hclust", combine_similarity_fun, number_of_clusters, distance_type)

# split cluster
clusters <- cutree(hclust_result[[3]], number_of_clusters)

