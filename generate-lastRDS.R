# generate-lastRDS.R
# raw IndividualWorak.RDS에서 출발해 필터링 파이프라인을 거쳐 lastRDS.rds 생성
#
# 파이프라인 순서:
#   길이 일관성 → 3D속도 → 공간 범위 → 관측 간격 → 개인정보 → 등산 시간

source("filtering-functions.R")

# ── 입력 경로 ─────────────────────────────────────────────────────────────────
TRAJ_RDS <- "C:/Users/erich/Dropbox/Worak/IndividualWorak.RDS"
INFO_CSV <- "C:/Users/erich/Dropbox/Worak/personalinfo.csv"
OUT_RDS  <- "C:/Users/erich/Desktop/학부/2024 학술제/lastRDS.rds"

# ── 파라미터 ──────────────────────────────────────────────────────────────────
SPEED_THRESHOLD    <- 10    # 최고 3D속도 상한 (km/h)
TIMEDIFF_THRESHOLD <- 1200  # GPS 관측 간격 상한 (초, 20분)
DURATION_THRESHOLD <- 12600 # 최소 등산 시간 (초, 3.5시간)

# ── 데이터 로드 ───────────────────────────────────────────────────────────────
original_RDS <- readRDS(TRAJ_RDS)
personalInfo <- read.csv(INFO_CSV, stringsAsFactors = FALSE)

# ── 파이프라인 ────────────────────────────────────────────────────────────────
cat(sprintf("원본 데이터: %d개\n", length(original_RDS)))

filtered1 <- filter_by_length_consistency(original_RDS)
cat(sprintf("길이 일관성 필터 후: %d개\n", length(filtered1)))

filtered2 <- filter_by_3Dspeed_threshold(filtered1, SPEED_THRESHOLD)
cat(sprintf("속도 필터 후 (<%dkm/h): %d개\n", SPEED_THRESHOLD, length(filtered2)))

filtered3 <- filter_by_spatial_boundary(filtered2)
cat(sprintf("공간 범위 필터 후: %d개\n", length(filtered3)))

filtered4 <- filter_by_timedifference(filtered3, TIMEDIFF_THRESHOLD)
cat(sprintf("관측 간격 필터 후 (<%d초): %d개\n", TIMEDIFF_THRESHOLD, length(filtered4)))

filtered5 <- filter_by_personal_information(filtered4, personalInfo)
cat(sprintf("개인정보 매칭 후: %d개\n", length(filtered5)))

filtered6 <- filter_by_duration_threshold(filtered5, DURATION_THRESHOLD)
cat(sprintf("등산 시간 필터 후 (>=%d초): %d개\n", DURATION_THRESHOLD, length(filtered6)))

# ── 저장 ──────────────────────────────────────────────────────────────────────
saveRDS(filtered6, OUT_RDS)
cat(sprintf("저장 완료: %s\n", OUT_RDS))
