library(ggplot2)
library(gridExtra)
library(tidyverse)
library(raster)
library(OpenStreetMap)
#library(hms) #time similarity 그림 그리는데 필요
library(lubridate)
library(sf)
library(hms)

RDS <- readRDS("C:/Users/erich/Dropbox/Worak/lastRDS.rds")

# 등산객들의 등산 바운더리 확인 ---------------------------------------------------------
find_map_boundary <- function(RDS){
  #등산객들의 등산 바운더리 확인
  eastest <- 0
  westest <- 10000
  northest <- 0
  southest <- 90
  for (i in seq(1,length(RDS))){
    east <- max(RDS[[i]]$longitude)
    west <- min(RDS[[i]]$longitude)
    north <- max(RDS[[i]]$latitude)
    south <- min(RDS[[i]]$latitude)
    
    if (east > eastest) {eastest <- east}
    if (west < westest) {westest <- west}
    if (north > northest) {northest <- north}
    if (south < southest) {southest <- south}
  }
  
  return(c(westest, southest, eastest, northest))
}

boundaries <- find_map_boundary(RDS)

# 위에서 확인한 등산 바운더리로 해당 부분 지도 가져오는 함수 ---------------------------------------------------------
find_map <- function(boundary_vector){
  #월악산 지도 확인
  map <- openmap(c(boundary_vector[4], boundary_vector[1]),
                 c(boundary_vector[2], boundary_vector[3]),
                 zoom = 10,
                 type = "bing",
                 mergeTiles = TRUE)
  map_p <- openproj(map)

  return(map_p)
}

map_p <- find_map(boundaries)

# 월악산 지도 그리는 함수 -------------------------------------------------------------------------
# 위에서 확인한 map_p를 통해,
# 월악산 지도를 그리는 함수.
# 투명도를 input으로 받아, 지도의 투명도를 결정정
draw_map <- function(alpha){
  p <- OpenStreetMap::autoplot.OpenStreetMap(map_p) +
    labs(x = "longitude", y = "latitude") +
    annotate(
      "rect",
      xmin = westest, xmax = eastest,   # 지도 범위에 맞춰 좌표 설정
      ymin = southest, ymax = northest,
      fill = "white", alpha = alpha    # 반투명 하얀색
    )
  
  return(p)
}

# 등산객의 등산 경로를 그려주는 함수 -------------------------------------------------------------------------
# indivList의 num번째 사람의 등산 경로를,
# alpha의 투명도를 가지는 color로 그림.
draw_route <- function(num, indivList, color, alpha){
  p <- geom_point(data = data.frame(indivList[[num]]),
                  aes(x = longitude, y = latitude), 
                  col = color, 
                  size = 1, 
                  alpha = alpha)
  
  return(p)
}