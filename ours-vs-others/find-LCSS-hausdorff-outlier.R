
# LCSS --------------------------------------------------------------------
library(geodist)

# 1) space_similarity: (주의) 첫 컬럼이 index면 제거
space_similarity <- read.csv("C:/Users/erich/Desktop/2024 학술제/repo/revision/space_sim_kdtree_chordal_eps0.1km.csv")
# space_similarity <- as.matrix(space_similarity[, -1])  # 필요 없으면 이 줄 제거

RDS <- readRDS("C:/Users/erich/Dropbox/Worak/lastRDS.rds")
N <- length(RDS)

# 2) start/end 좌표 행렬 만들기 (N x 2: lon, lat)
starts <- cbind(
  lon = sapply(RDS, function(x) x$longitude[1]),
  lat = sapply(RDS, function(x) x$latitude[1])
)

ends <- cbind(
  lon = sapply(RDS, function(x) tail(x$longitude, 1)),
  lat = sapply(RDS, function(x) tail(x$latitude, 1))
)

# 3) 거리 행렬 3개만 계산하면 됨
# D_ss: start-start
D_ss <- geodist(starts, starts, measure = "geodesic")  # N x N
# D_ee: end-end
D_ee <- geodist(ends, ends, measure = "geodesic")      # N x N
# D_se: start-end
D_se <- geodist(starts, ends, measure = "geodesic")    # N x N

# 4) direct / reverse (meters)
direct  <- D_ss + D_ee
reverse <- D_se + t(D_se)   # end-start는 start-end의 전치

identical(D_se, t(D_es)) # 확인용

ratio <- direct / (reverse + 1e-9)

# 5) upper triangle만 뽑아서 결과 테이블 만들기
ut <- upper.tri(ratio, diag = FALSE)
idx <- which(ut, arr.ind = TRUE)

result <- data.frame(
  i = idx[,1],
  j = idx[,2],
  space_similarity = space_similarity[ut],
  direct_m = direct[ut],
  reverse_m = reverse[ut],
  ratio = ratio[ut]
)

# 6) 원하는 기준으로 정렬해서 “LCSS 낮을 것 같은 후보” 상위 확인
# (반대방향 강할수록 ratio가 커짐)
result_sorted <- result[order(-result$ratio), ]
head(result_sorted, 50)

write.csv(result_sorted, "C:/Users/erich/Desktop/2024 학술제/similarity-calculation/space-similarity-ratio.csv", row.names = FALSE)


# Hausdorff ---------------------------------------------------------------
hausdorff <- read.csv("C:/Users/erich/Desktop/2024 학술제/similarity-calculation/hausdorff_dist_utm.csv")

# result_sorted에 hausdorff_measure 추가
result_sorted$hausdorff_measure <- hausdorff[cbind(result_sorted$i, result_sorted$j)]

summary(result_sorted$hausdorff_measure)

plot(result_sorted$space_similarity,
     result_sorted$hausdorff_measure,
     pch = 16, cex = 0.4,
     xlab = "Space Similarity",
     ylab = "Hausdorff Distance (m)")

abline(v = 0.95, col = "red")
abline(h = h_cut, col = "blue")


subset(result_sorted,
       space_similarity < 0.1 &
         hausdorff_measure < 2000)


