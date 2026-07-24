# Worak
Partial source codes for Worak data analysis

- `20250215 Similarities.R`: Similarity matrix calculation for the data analysis

- `20250216 Clustering.R`: Clustering analysis based on the similarity matrix

- Similarity outputs are `space_sim_939.csv`, `time_sim_939.csv`, and `velocity_sim_939.csv`

## Generating lastRDS.rds

`lastRDS.rds` is the filtered dataset used for similarity and clustering analysis.

**Required input files:**
- `IndividualWorak.RDS` — raw trajectory data
- `personalinfo.csv` — hikers' personal information

**Steps:**

1. Open `generate-lastRDS.R` and set the three path variables at the top:
   ```r
   TRAJ_RDS <- "write/your/path/to/IndividualWorak.RDS"
   INFO_CSV <- "write/your/path/to/personalinfo.csv"
   OUT_RDS  <- "write/your/path/to/lastRDS.rds"
   ```

2. Run from the repo directory:
   ```bash
   Rscript generate-lastRDS.R
   ```

**Filtering pipeline** (applied in order):
1. Length consistency — removes trajectories with mismatched field lengths
2. 3D speed < 10 km/h — removes outlier speeds
3. Spatial boundary — keeps trajectories within Woraksan area
4. Observation interval < 1200 s — removes trajectories with large GPS gaps
5. Personal info match — keeps hikers with valid height/weight/gender records
6. Duration ≥ 12600 s — keeps hikes of at least 3.5 hours