library(ggplot2)
library(gridExtra)
library(tidyverse)
library(raster)
library(OpenStreetMap)
#library(hms) #time similarity 그림 그리는데 필요
library(lubridate)
library(sf)
library(hms)
library(patchwork)
library(scales)
library(patchwork)

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

# 두 trajectory의 유사도 시각화 함수 -------------------------------------------------------------------------
# 유사도 시각화 -----------------------------------------------------------------
# indivList의 num1, num2의 유사도를 시각화 하는 함수입니다.
similarity_visualize <- function(num1, 
                                 num2, 
                                 RDS, 
                                 b_mark_start_and_end = TRUE) {
  
  # figure settings ---------------------------------------------------------
  windowsFonts(TimesNewer = windowsFont("Times Newer Roman"))
  
  big_theme <- theme_minimal(base_size = 16) +
    theme(
      text = element_text(family = "TimesNewer"),
      axis.title = element_text(size = 18, face = "bold", family = "TimesNewer"),
      axis.text  = element_text(size = 15, family = "TimesNewer"),
      legend.title = element_text(size = 16, family = "TimesNewer"),
      legend.text  = element_text(size = 14, family = "TimesNewer")
    )
  
  # (a) Spatial ----------------------------------------------------------------
  spatial <- OpenStreetMap::autoplot.OpenStreetMap(map_p) +
    geom_point(data = data.frame(RDS[[num1]]),
               aes(x = longitude, y = latitude), col = "red", alpha = 0.5) +
    geom_point(data = data.frame(RDS[[num2]]),
               aes(x = longitude, y = latitude), col = "blue", alpha = 0.5) +
    labs(caption = "(a) Space Similarity",
         x = "longitude",
         y = "latitude") +
    scale_x_continuous(n.breaks = 2) + 
    big_theme +
    theme(
      plot.caption = element_text(size = 16, family = "TimesNewer", hjust = 0.5)
    )
  
  # (b) Spatial Detailed ---------------------------------------------------------
  traj1 <- data.frame(RDS[[num1]])
  traj2 <- data.frame(RDS[[num2]])
  
  detailed_spatial <- ggplot() +
    geom_point(data = traj1, aes(x = longitude, y = latitude), col = "red", alpha = 0.5) +
    geom_point(data = traj2, aes(x = longitude, y = latitude), col = "blue", alpha = 0.5) +
    labs(caption = "(b) Space Similarity (Detail)",
         x = "longitude",
         y = "latitude") +
    scale_x_continuous(n.breaks = 2) +
    big_theme +
    theme(
      plot.caption = element_text(size = 16, family = "TimesNewer", hjust = 1)
    )
  
  # add start/end markers only if requested
  if (isTRUE(b_mark_start_and_end)) {
    
    # offset 계산 (데이터 범위의 2%)
    x_range <- range(c(traj1$longitude, traj2$longitude), na.rm = TRUE)
    y_range <- range(c(traj1$latitude,  traj2$latitude),  na.rm = TRUE)
    
    x_offset <- diff(x_range) * 0.01
    y_offset <- diff(y_range) * 0.01
    
    se1 <- rbind(
      data.frame(longitude = traj1$longitude[1],           latitude = traj1$latitude[1],           label = "Start"),
      data.frame(longitude = traj1$longitude[nrow(traj1)], latitude = traj1$latitude[nrow(traj1)], label = "End")
    )
    se2 <- rbind(
      data.frame(longitude = traj2$longitude[1],           latitude = traj2$latitude[1],           label = "Start"),
      data.frame(longitude = traj2$longitude[nrow(traj2)], latitude = traj2$latitude[nrow(traj2)], label = "End")
    )
    
    detailed_spatial <- detailed_spatial +
      # 빨강: 텍스트 오른쪽 위
      geom_point(data = se1, aes(x = longitude, y = latitude),
                 col = "red", size = 3, alpha = 0.9) +
      geom_text(data = se1, aes(x = longitude, y = latitude, label = label),
                color = "red",
                nudge_x = x_offset,
                nudge_y = y_offset,   # 오른쪽 위
                hjust = 0, vjust = 0,
                size = 5, fontface = "bold") +
      
      # 파랑: 텍스트 오른쪽 아래
      geom_point(data = se2, aes(x = longitude, y = latitude),
                 col = "blue", size = 3, alpha = 0.9) +
      geom_text(data = se2, aes(x = longitude, y = latitude, label = label),
                color = "blue",
                nudge_x = x_offset,
                nudge_y = -y_offset,  # 오른쪽 아래
                hjust = 0, vjust = 1,
                size = 5, fontface = "bold")
  }
  
  
  # Velocity ----------------------------------------------------------------
  df1 <- data.frame(value = indivList[[num1]]$spd3D[-1], group = as.character(num1))
  df2 <- data.frame(value = indivList[[num2]]$spd3D[-1], group = as.character(num2))
  df <- rbind(df1, df2)
  
  velocity <- ggplot(df, aes(x = value, y = ..density.., fill = group)) +
    geom_density(alpha = 0.5) + 
    scale_fill_manual(values = c("red", "blue")) +
    scale_x_continuous(limits = c(0, 10)) +
    labs(caption = "(c) Velocity Overlapping",
         x = "spd3D (km/h)",
         y = "Frequency") +
    big_theme +
    theme(
      legend.position = "none",
      plot.caption = element_text(size = 16, family = "TimesNewer", hjust = 0.5)
    )
  
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
    scale_x_time(
      breaks = scales::breaks_width("1 hour"), 
      labels = scales::label_time(format = "%H")
    ) +
    labs(
      caption = "(d) Time Overlapping",
      x = "Time (hour)",
      y = "Group"
    ) +
    coord_cartesian(ylim = c(-10, 10)) +
    scale_y_continuous(labels = NULL) +
    big_theme +
    theme(
      legend.position = "none",
      plot.caption = element_text(size = 16, family = "TimesNewer", hjust = 0.5)
    )
  
  plot_list <- list(
    spatial = spatial,
    detailed_spatial = detailed_spatial,
    velocity = velocity,
    time = time
  )
  
  return(plot_list)
}



# 월악산 지도 그리는 함수 -------------------------------------------------------------------------
# 위에서 확인한 map_p를 통해,
# 월악산 지도를 그리는 함수.
# 투명도를 input으로 받아, 지도의 투명도를 결정정
draw_map <- function(alpha, map_p){
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

boundaries <- find_map_boundary(RDS)
map_p <- find_map(boundaries)

plots <- similarity_visualize(90, 1, RDS)

# png("90 vs 1 new.png", width = 600, height = 630)
grid.arrange(
  plots$spatial, plots$detailed_spatial,
  nullGrob(), nullGrob(),  # 🔹 빈 행
  plots$velocity, plots$time,
  nrow = 3, ncol = 2,
  heights = c(1, 0.05, 1)  # 🔹 가운데 행 높이 조절
)
# dev.off()
