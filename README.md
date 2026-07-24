# Worak
월악산 등산 데이터 분석 소스 코드

- `20250215 Similarities.R`: 유사도 행렬 계산
- `20250216 Clustering.R`: 유사도 행렬 기반 클러스터링 분석
- 유사도 결과물: `space_sim_939.csv`, `time_sim_939.csv`, `velocity_sim_939.csv`

## lastRDS.rds 생성

`lastRDS.rds`는 유사도 계산 및 클러스터링에 사용되는 필터링된 데이터셋입니다.

**필요한 입력 파일:**
- `IndividualWorak.RDS` — 등산객 궤적 원본 데이터
- `personalinfo.csv` — 등산객 개인정보

**실행 방법:**

1. `generate-lastRDS.R` 상단의 경로 변수 3개를 설정합니다:
   ```r
   TRAJ_RDS <- "write/your/path/to/IndividualWorak.RDS"
   INFO_CSV <- "write/your/path/to/personalinfo.csv"
   OUT_RDS  <- "write/your/path/to/lastRDS.rds"
   ```

2. repo 디렉토리에서 실행합니다:
   ```bash
   Rscript generate-lastRDS.R
   ```

**필터링 파이프라인** (순서대로 적용):
1. 길이 일관성 — longitude/latitude/elev/tm 관측 수가 다른 궤적 제거
2. 3D 속도 < 10 km/h — 속도 이상치 제거
3. 공간 범위 — 월악산 영역 밖 궤적 제거
4. 관측 간격 < 1200초 — GPS 시간 간격이 큰 궤적 제거
5. 개인정보 매칭 — 키/몸무게/성별 정보가 유효한 등산객만 유지
6. 등산 시간 ≥ 12600초 — 3.5시간 미만 등산 제거
