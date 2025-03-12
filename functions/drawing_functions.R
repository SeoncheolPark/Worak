library(ggplot2)
library(gridExtra)
library(tidyverse)
library(raster)
library(OpenStreetMap)
#library(hms) #time similarity 그림 그리는데 필요
library(lubridate)
library(sf)
library(hms)

df <- readRDS("path/to/your/rds.rds")

# 등산객들의 등산 바운더리 확인 ---------------------------------------------------------
find_map_boundary <- function(indivList){
  #등산객들의 등산 바운더리 확인
  eastest <- 0
  westest <- 10000
  northest <- 0
  southest <- 90
  for (i in seq(1,length(indivList))){
    east <- max(indivList[[i]]$longitude)
    west <- min(indivList[[i]]$longitude)
    north <- max(indivList[[i]]$latitude)
    south <- min(indivList[[i]]$latitude)
    
    if (east > eastest) {eastest <- east}
    if (west < westest) {westest <- west}
    if (north > northest) {northest <- north}
    if (south < southest) {southest <- south}
  }
  
  return(c(westest, southest, eastest, northest))
}

# 위에서 확인한 등산 바운더리로 해당 부분 지도 가져오는 함수 ---------------------------------------------------------
find_map <- function(boundary_vector){
  #월악산 지도 확인
  map <- openmap(c(boundary_vector[4], boundary_vector[1]), 
                 c(boundary_vector[2], boundary_vector[3]), 
                 zoom = 10, 
                 type = "osm", 
                 mergeTiles = TRUE)
  map_p <- openproj(map)
  
  return(map_p)
}

map_p <- find_map(find_map_boundary(df))

# 두 trajectory의 유사도 시각화 함수 -------------------------------------------------------------------------
# 유사도 시각화 -----------------------------------------------------------------
similarity_visualize <- function(num1, num2, indivList) {
  # Spatial ----------------------------------------------------------------
  spatial <- OpenStreetMap::autoplot.OpenStreetMap(map_p) +
    geom_point(data = data.frame(indivList[[num1]]),
               aes(x = longitude, y = latitude), col = "red", alpha = 0.5) +
    geom_point(data = data.frame(indivList[[num2]]),
               aes(x = longitude, y = latitude), col = "blue", alpha = 0.5) +
    labs(caption = paste("(a) Spatial Similarity of",num1,"vs",num2, sep =" "), 
         x = "longitude",
         y = "latitude") +
    theme(plot.caption = element_text(hjust = 0.5, size = 20))
  
  # Spatial Detailed ---------------------------------------------------------
  detailed_spatial <- ggplot() +
    geom_point(data = data.frame(indivList[[num1]]),
               aes(x = longitude, y = latitude), col = "red", alpha = 0.5) +
    geom_point(data = data.frame(indivList[[num2]]),
               aes(x = longitude, y = latitude), col = "blue", alpha = 0.5) +
    theme_minimal() + 
    labs(caption = paste("(b) Spatial Similarity of",num1,"vs",num2, sep =" "), 
         x = "longitude",
         y = "latitude") +
    theme(plot.caption = element_text(hjust = 0.5, size = 20))
  
  # Velocity ----------------------------------------------------------------
  df1 <- data.frame(value = indivList[[num1]]$spd3D[-1], group = as.character(num1))
  df2 <- data.frame(value = indivList[[num2]]$spd3D[-1], group = as.character(num2))
  df <- rbind(df1, df2)
  
  # velocity <- ggplot(df, aes(x = value, y = ..density.., fill = group)) +
  #   geom_histogram(position = "identity", alpha = 0.5, binwidth = 1) +
  #   scale_fill_manual(values = c("red", "blue")) +
  #   labs(caption = paste("(c) Velocity Overlapping of",num1,"vs",num2, sep =" "),
  #        x = "spd3D",
  #        y = "Frequency") +
  #   theme_minimal() +
  #   theme(plot.caption = element_text(hjust = 0.5, size = 20))

  velocity <- ggplot(df, aes(x = value, y = ..density.., fill = group)) +
    # geom_histogram(position = "identity", alpha = 0.5, binwidth = 1) +
    geom_density(alpha = 0.5) + 
    scale_fill_manual(values = c("red", "blue")) +
    scale_x_continuous(limits = c(0, 10)) +
    labs(caption = paste("(c) Velocity Overlapping of",num1,"vs",num2, sep =" "),
         x = "spd3D",
         y = "Frequency") +
    theme_minimal() +
    theme(plot.caption = element_text(hjust = 0.5, size = 20))

  # 
  # Time ----------------------------------------------------------------
  time1 <- format(indivList[[num1]]$time, '%H:%M:%S')
  time1 <- as_hms(time1)
  df1 <- data.frame(time = time1, y = 1, group = as.character(num1))
  time2 <- format(indivList[[num2]]$time, '%H:%M:%S')
  time2 <- as_hms(time2)
  df2 <- data.frame(time = time2, y = 2, group = as.character(num2))
  
  df <- rbind(df1, df2)
  
  time <- ggplot(df, aes(x = time, y = y, color = group)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = c("red", "blue")) +
    labs(caption = paste("(d) Time Overlapping of",num1,"vs",num2, sep =" "),
         x = "Time",
         y = "Group") +
    theme_minimal() +
    coord_cartesian(ylim = c(-10,10)) + 
    scale_y_continuous(labels = NULL) + 
    theme(axis.text.y = element_text(size = 12),  #y축 글자크기 조정
          plot.caption = element_text(hjust = 0.5, size = 20)) 
  
  # p <- grid.arrange(spatial, velocity, time, nrow = 1, ncol = 3)
  # return(p)
  plot_list <- list(spatial = spatial,
                    detailed_spatial = detailed_spatial,
                    velocity = velocity,
                    time = time)
  
  return(plot_list)
}


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

# 한 trajectory의 값들 시각화 -------------------------------------------------------------------------
# 한 trajectory의 long, lat, elev, spd 등을 시각화 하는 함수 입니다.
# x축과 y축으로 쓰일 값을 입력 받아, scatter plot을 그려줍니다.

draw_trajectory_eda <- function(indiv_df, col1, col2, fig_name){
  indiv_df <- data.frame(indiv_df)
  indiv_df <- data.frame(x = indiv_df[[col1]], y = indiv_df[[col2]])
  
  p <- ggplot(data = indiv_df, aes(x = x, y = y)) + 
    geom_point() +
    theme_minimal() +
    labs(x = col1, y = col2, caption = fig_name) +
    theme(plot.cpation = element_text(hjust = 0.5, size = 20),
          axis.title.x = element_text(size = 15),   # x축 제목 폰트 크기
          axis.title.y = element_text(size = 15),   # y축 제목 폰트 크기
          axis.text.x = element_text(size = 12),    # x축 눈금 폰트 크기
          axis.text.y = element_text(size = 12)     # y축 눈금 폰트 크기
    ) 
  
  return(p)
}

indivNum <- 12 # 그림 그릴 개인의 번호
spd3D <- draw_trajectory_eda(df[[indivNum]], "time", "spd3D", "(f) spd3D")
spd2D <- draw_trajectory_eda(df[[indivNum]], "time", "spd2D", "(e) spd2D")
latitude <- draw_trajectory_eda(df[[indivNum]], "time", "latitude", "(b) Latitude")
longitude <- draw_trajectory_eda(df[[indivNum]], "time", "longitude", "(a) Longitude")
elev <- draw_trajectory_eda(df[[indivNum]], "time", "elev", "(c) Elevation")
longlat <- draw_trajectory_eda(df[[indivNum]], "longitude", "latitude", "(d) Long. vs Lat.")

png("figure_file_name.png", width = 1200, height = 800)
grid.arrange(
  longitude, latitude, elev,
  longlat, spd2D, spd3D,
  nrow = 2, ncol = 3 # 2행, 3열로 배치
)
dev.off()