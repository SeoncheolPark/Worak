# Similarity --------------------------------------------------------------

## 1. Space -------------------------------------------------------------------
library(geodist)
similarity_spatial <- function(obj1, obj2, epsilon=0.1){
  #obj1 and obj2 contain spatial trajectories
  #epsilon: tuning parameter, 한 point를 기준으로 반지름이 epsilon인 ball을
  #생각하고 그 안에 다른 obj의 point가 들어있는지 없는지를 체크한다.
  
  #calculate lengths
  n1 <- length(obj1$long)
  n2 <- length(obj2$long)
  
  LongLat_matrix_obj1 <- data.frame(lon=obj1$long, lat=obj1$lat)
  LongLat_matrix_obj2 <- data.frame(lon=obj2$long, lat=obj2$lat)
  
  #학생들이 만들어 놓은 chordal_dist 함수 이용
  chordal_dist_vec=function(x,y){
    #01. Calculate geodesic distance first
    geodesic_matrix <- geodist(x,y, measure = "geodesic")/1000 #거리 km단위로 변환
    
    #02. Convert geodesic distance into chordal distance(Lecture01 p16, p17)
    radius <- 6371 #지구 반지름의 평균값(단위: km)
    theta <- geodesic_matrix/radius
    chordal_matrix <- 2 * radius * sin(theta/2)
    
    return(chordal_matrix) #diagonal에서 딱 하나 넘어간게 두 시점간의 거리(단위 km)
  }
  #(1) obj1에 대해 계산
  #result_obj1: n2*n1 matrix
  result_obj1 <- apply(LongLat_matrix_obj1, 1, function(x) chordal_dist_vec(x, LongLat_matrix_obj2))
  is_neighbor_obj1 <- apply(result_obj1, 2, function(x) sum(x<epsilon))
  result_vec1 <- is_neighbor_obj1>0
  result1 <- sum(result_vec1)/length(result_vec1)
  
  #(2) obj2에 대해 계산
  #n1 < n2인 경우, obj2의 모든 gps observation 중 obj1에 가까이 있는 점의 갯수를 센다
  result_obj2 <- apply(LongLat_matrix_obj2, 1, function(x) chordal_dist_vec(x, LongLat_matrix_obj1))
  is_neighbor_obj2 <- apply(result_obj2, 2, function(x) sum(x<epsilon))
  result_vec2 <- is_neighbor_obj2>0
  result2 <- sum(result_vec2)/length(result_vec2)
  
  similarity_list=list(obj1=result_vec1, obj2=result_vec2)
  
  #두 개의 result중 작은 값을 similarity measure로 삼는다)
  return(list(similarity=min(result1, result2), similarity_list=similarity_list, similarity_obj1=result1, similarity_obj2=result2))
}


## 2. Velocity -------------------------------------------------------------
library(overlapping)
similarity_velocity <- function(obj1, obj2, plot = FALSE){
  X <- list(X1 = obj1$spd3D[-1],
            X2 = obj2$spd3D[-1])
  velocity_similarity <- overlap(X, type="1", plot = plot)$OV
  return(velocity_similarity)
}

## 3. Time -----------------------------------------------------------------
#install.packages("lubridate")
library(lubridate)

#obj1 <- individualWorak[[i]]
#obj2 <- individualWorak[[j]]
convert_to_seconds <- function(HOUR, MIN, SECOND){
  #시간을 초(numeric)로 바꾸어주는 함수
  return(HOUR * 60 * 60 + MIN * 60 + SECOND)
}

convert_to_time <- function(secondsInNumeric){
  #초(numeric)를 시간으로 바꾸어 주는 함수
  HOUR <- secondsInNumeric %/% 3600
  secondsInNumeric <- secondsInNumeric %% 3600
  MIN <- secondsInNumeric %/% 60
  secondsInNumeric <- secondsInNumeric %% 60
  SECOND <- secondsInNumeric
  return(paste(HOUR,MIN,SECOND, sep=":"))
}

similarity_time <- function(obj1, obj2){
  #두 obj의 time_similarity 계산
  
  #계산을 위한 사전 준비
  #계산의 IDEA
  # 1) 어차피 날짜는 큰 상관이 없음
  # 2) 시간 자료에서 시간:분:초 만 추출
  # 3) obj1의 등산 시작 시간 ~ 종료 시간과 obj2의 등산 시작 ~ 종료 시간의
  #    겹치는 넓이 확인
  # +) 이때 obj1과 obj2의 총 등산 시간에 따라, 겹치는 부분의 비율이 달라짐
  #    그래서 similarity_obj1 과 similarity_obj2를 따로 계산
  #    spatial similarity를 계산할 때 최솟값을 return 한 것과 같이
  #    여기서도 최솟값을 return
  
  obj1_start <- obj1$time[1]
  obj1_end <- obj1$time[length(obj1$time)]
  
  adjust <- 0
  if (day(obj1_start) != day(obj1_end)){
    adjust <- 24
  }
  
  obj1_start <- convert_to_seconds(hour(obj1_start), minute(obj1_start), second(obj1_start))
  obj1_end <- convert_to_seconds(hour(obj1_end)+adjust, minute(obj1_end), second(obj1_end))
  
  obj1_times <- seq(obj1_start, obj1_end)
  
  obj2_start <- obj2$time[1]
  obj2_end <- obj2$time[length(obj2$time)]
  
  adjust <- 0
  if (day(obj2_start) != day(obj2_end)){
    adjust <- 24
  }
  
  obj2_start <- convert_to_seconds(hour(obj2_start), minute(obj2_start), second(obj2_start))
  obj2_end <- convert_to_seconds(hour(obj2_end)+adjust, minute(obj2_end), second(obj2_end))
  
  obj2_times <- seq(obj2_start, obj2_end)
  
  timeOverlap <- intersect(obj1_times, obj2_times)
  
  if (length(timeOverlap) > 0 ){ #겹치는 시간이 있는 경우만 체크
    #겹치기 시작한 시간과, 겹치는게 끝난 시간
    overlapStart <- convert_to_time(timeOverlap[1])
    overlapEnd <- convert_to_time(timeOverlap[length(timeOverlap)])
    
    timeOverlap <- timeOverlap[length(timeOverlap)] - timeOverlap[1]
    #겹치는 시간이 얼마나 길었는지?
    howlong <- convert_to_time(timeOverlap)
    
    similarity_time_1 <- timeOverlap/length(obj1_times)
    similarity_time_2 <- timeOverlap/length(obj2_times)
    return(list(similarity = min(similarity_time_1,similarity_time_2),
                similarity_obj1 = similarity_time_1,
                similarity_obj2 = similarity_time_2,
                timeOverlap = c(overlapStart, overlapEnd, howlong)))
  }
  
  else{
    return(list(similarity = 0,
                similarity_obj1 = 0,
                similarity_obj2 = 0,
                timeOverlap = c("0", "0", "0")))
  }
}


# Calculate Similarity ----------------------------------------------------
library(foreach)
library(doParallel)

df <- readRDS("path/to/lastRDS.rds")
data_length <- len(df) #939

space_sim <- data.frame(matrix(0, nrow = data_length, ncol = data_length))
velocity_sim <- data.frame(matrix(0, nrow = data_length, ncol = data_length))
time_sim <- data.frame(matrix(0, nrow = data_length, ncol = data_length))

#병렬 처리 준비
numCores <- detectCores()/2           # 시스템의 전체 코어 수 - 1
cl <- parallel::makeCluster(numCores)    # 클러스터 생성
doParallel::registerDoParallel(cl)       # 병렬 백엔드 등록



## Space -------------------------------------------------------------------
space_sim <- foreach(i = 1:data_length, 
                       .packages = c("geodist"),
                       .combine = rbind) %dopar% {
                         foreach(j = 1:data_length,
                                 .combine = "c") %do% {
                                   p1 <- df[[i]]
                                   p2 <- df[[j]]
                                   
                                   similarity_spatial(p1, p2)$similarity_obj1
                                 }
                       }

write.csv(space_sim, "path/to/save/space_sim.csv")


## Velocity ----------------------------------------------------------------
velocity_sim <- foreach(i = 1:data_length, 
                        .packages = c("foreach","overlapping"),
                        .combine = rbind) %dopar% {
                          foreach(j = 1:data_length,
                                  .combine = "c") %do% {
                                    p1 <- df[[i]]
                                    p2 <- df[[j]]
                                    
                                    similarity_velocity(p1, p2)
                                  }
                        }

velocity_sim <- data.frame(velocity_sim)
write.csv(velocity_sim, "path/to/save/velocity_sim.csv")


## Time --------------------------------------------------------------------
time_sim <- foreach(i = 1:data_length,
                    .packages = c("foreach","lubridate"),
                    .combine = rbind) %dopar% {
                      foreach(j = 1:data_length,
                              .combine = "c") %do% {
                                p1 <- imsi[[i]]
                                p2 <- imsi[[j]]
                                
                                similarity_time(p1, p2)$similarity_obj1
                              }
                    }

write.csv(time_sim, "path/to/save/time_sim.csv")