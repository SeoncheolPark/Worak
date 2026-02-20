# Ours New vs LCSS ------------------------------------------------------
# drawing/basic-functions.R 함수 실행 필요
boundaries <- find_map_boundary(RDS)
map_p <- find_map(boundaries)


## Ours와 LCSS가 큰 차이 나는 경우 --------------------------------------------------
# similarity-calculation/space-similarity-ratio.csv 파일 확인.
ours_vs_LCSS <- function(A, B){
  plots <- similarity_visualize(A, B, RDS)
  plots$detailed_spatial
}

A <- 276
B <- 749

ours_vs_LCSS(A, B)

# Ours New vs Hausdorff ---------------------------------------------------
hausdorff <- read.csv("C:/Users/erich/Desktop/2024 학술제/similarity-calculation/hausdorff_dist_utm.csv")

ours_vs_hausdorff <- function(A, B){
  plots <- similarity_visualize(A, B, RDS)
  grid.arrange(
    plots$spatial, plots$detailed_spatial,
    nrow = 1, ncol = 2
  )
}

A <- 561
B <- 588

ours_vs_hausdorff(A, B)