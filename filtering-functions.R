library(dplyr)

filter_by_length_consistency <- function(RDS,
                                         required_fields = c("longitude","latitude","elev","tm")) {
  keep <- logical(length(RDS))

  for (person_id in seq_along(RDS)) {
    traj <- RDS[[person_id]]

    if (is.null(traj)) { keep[person_id] <- FALSE; next }
    if (!all(required_fields %in% names(traj))) { keep[person_id] <- FALSE; next }

    lengths <- c(length(traj$longitude), length(traj$latitude),
                 length(traj$elev),      length(traj$tm))

    if (all(lengths == 0L) || any(lengths == 0L)) { keep[person_id] <- FALSE; next }
    if (length(unique(lengths)) != 1L)             { keep[person_id] <- FALSE; next }

    keep[person_id] <- TRUE
  }

  return(RDS[keep])
}

filter_by_3Dspeed_threshold <- function(RDS, threshold) {
  result <- list()
  idx <- 1

  for (i in seq_along(RDS)) {
    spd3Dvec <- as.vector(RDS[[i]]$spd3D[-1])
    if (anyNA(spd3Dvec) || any(is.nan(spd3Dvec))) next
    if (max(spd3Dvec) > threshold) next
    result[[idx]] <- RDS[[i]]
    idx <- idx + 1
  }

  return(result)
}

filter_by_spatial_boundary <- function(
    RDS,
    northest = 37,
    eastest  = 128.4,
    southest = 36.5,
    westest  = 128
) {
  keep_index <- sapply(RDS, function(traj) {
    if (is.null(traj)) return(FALSE)
    if (!all(c("longitude","latitude") %in% names(traj))) return(FALSE)
    if (anyNA(traj$longitude) || anyNA(traj$latitude)) return(FALSE)
    all(traj$latitude  > southest & traj$latitude  < northest &
        traj$longitude > westest  & traj$longitude < eastest)
  })

  return(RDS[keep_index])
}

filter_by_timedifference <- function(RDS, threshold.second) {
  keep_index <- sapply(RDS, function(traj) {
    if (is.null(traj)) return(FALSE)
    if (!("tmdiff" %in% names(traj))) return(FALSE)
    if (anyNA(traj$tmdiff)) return(FALSE)
    all(traj$tmdiff < threshold.second)
  })

  return(RDS[keep_index])
}

filter_by_personal_information <- function(RDS, personalInfo) {
  stopifnot(is.list(RDS), is.data.frame(personalInfo))

  required_cols <- c("트랙명", "회원키값", "회원체중값", "회원성별구분코드")
  missing_cols <- setdiff(required_cols, colnames(personalInfo))
  if (length(missing_cols) > 0)
    stop("personalInfo에 필요한 컬럼이 없습니다: ", paste(missing_cols, collapse = ", "))

  track_names <- as.character(personalInfo$트랙명)
  n <- nchar(track_names)
  personalInfo$privateNum <- ifelse(
    is.na(track_names) | n < 10,
    NA_character_,
    substr(track_names, 6, n - 4)
  )

  personalInfo_f <- subset(
    personalInfo,
    !(회원키값 == 0 | 회원체중값 == 0) &
      (회원성별구분코드 %in% c("남", "여")) &
      (회원키값 >= 140 & 회원키값 < 200) &
      !is.na(privateNum)
  )

  validNames <- unique(as.character(personalInfo_f$privateNum))

  rds_names <- vapply(RDS, function(x) {
    if (is.null(x) || is.null(x$name)) NA_character_ else as.character(x$name)
  }, FUN.VALUE = character(1))

  return(RDS[!is.na(rds_names) & (rds_names %in% validNames)])
}

filter_by_duration_threshold <- function(RDS, DTT) {
  result <- list()
  idx <- 1

  for (i in seq_along(RDS)) {
    totalDT <- RDS[[i]]$tm[length(RDS[[i]]$tm)]
    if (totalDT >= DTT) {
      result[[idx]] <- RDS[[i]]
      idx <- idx + 1
    }
  }

  return(result)
}
