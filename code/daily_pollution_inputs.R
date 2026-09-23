# Published daily inputs for the revised baseline and time-varying analyses.
daily_pollution_releases <- data.frame(
  pollutant = c("pm25", "o3"),
  tag = c("lghap-pm25-zcta-daily-v1", "o3-zcta-daily-v1"),
  directory = c("lghap_pm25_zcta_daily_parquet", "o3_zcta_daily_parquet"),
  prefix = c("lghap_pm25_zcta_daily_", "o3_zcta_daily_"),
  value_column = c("pm25_ug_m3", "o3_ppb")
)
daily_pollution_end_date <- as.Date("2024-12-31")

trim_monthly_cohort <- function(cohort) {
  # Avoid replicating unused baseline labs and exposure audit columns every month.
  columns <- c("PX_ID", "PERS_ID", "waitlist_row_id", "WL_ORG", "candidate_zip",
               "index_date", "observed_end_date", "adverse_event", "age", "sex", "race",
               "zcta_svi_proxy", "listing_center", "listing_year", "organ_score",
               "kidney_dialysis_years", "kidney_no_dialysis_time", "kidney_diabetes",
               "age_group", "race_group", "multi_organ_candidate")
  data.table::copy(cohort[, intersect(columns, names(cohort)), with = FALSE])
}

expand_monthly_cohort <- function(intervals) {
  # Calendar lookup avoids constructing a separate Date sequence for every candidate.
  start_id <- as.integer(format(intervals$index_date, "%Y")) * 12L +
    as.integer(format(intervals$index_date, "%m")) - 1L
  end_id <- as.integer(format(intervals$analysis_end_date, "%Y")) * 12L +
    as.integer(format(intervals$analysis_end_date, "%m")) - 1L
  lengths <- end_id - start_id + 1L
  ids <- rep(start_id, lengths) + sequence(lengths) - 1L
  calendar_ids <- seq.int(min(start_id), max(end_id) + 1L)
  calendar <- as.Date(sprintf("%04d-%02d-01", calendar_ids %/% 12L, calendar_ids %% 12L + 1L))
  out <- intervals[rep(seq_len(nrow(intervals)), lengths)]
  out[, month_start := calendar[ids - min(start_id) + 1L]]
  out[, month_end := calendar[ids - min(start_id) + 2L] - 1L]
  out
}

daily_pollution_files <- function(path) {
  spec <- daily_pollution_releases[
    daily_pollution_releases$directory == basename(path), , drop = FALSE
  ]
  if (nrow(spec) != 1L) stop("Unrecognized daily pollution directory: ", path)
  files <- file.path(path, paste0(spec$prefix, 2005:2024, ".parquet"))
  if (any(!file.exists(files))) {
    stop("Missing published daily inputs: ", paste(basename(files[!file.exists(files)]), collapse = ", "),
         ". Obtain from https://github.com/petergraffy/environment_transplant_survival/releases/tag/", spec$tag)
  }
  files
}

daily_pollution_cache_path <- function(cache_file, files, value_col, out_col, resolution) {
  info <- file.info(files)
  signature <- tempfile()
  on.exit(unlink(signature))
  saveRDS(list(version = 1L, files = normalizePath(files, winslash = "/"),
               size = info$size, modified = as.numeric(info$mtime),
               value = value_col, output = out_col, resolution = resolution), signature)
  key <- unname(tools::md5sum(signature))
  sub("[.]csv[.]gz$", paste0("_", key, ".csv.gz"), cache_file)
}

read_daily_pollution_aggregate <- function(path, value_col, out_col,
                                         resolution = c("annual", "monthly"),
                                         cache_file = NULL) {
  resolution <- match.arg(resolution)
  files <- daily_pollution_files(path)
  if (is.null(cache_file)) {
    cache_file <- file.path("data", "cache", "daily_pollution",
                            paste0(value_col, "_", resolution, ".csv.gz"))
  }
  cache_file <- daily_pollution_cache_path(cache_file, files, value_col, out_col, resolution)
  if (file.exists(cache_file)) {
    return(readr::read_csv(cache_file, col_types = readr::cols(zip = readr::col_character()),
                           show_col_types = FALSE))
  }
  # Process one year at a time to bound memory use with national daily data.
  groups <- if (resolution == "monthly") c("zip", "year", "month") else c("zip", "year")
  minimum_days <- if (resolution == "monthly") 20L else 300L
  pieces <- lapply(files, function(file) {
    message("Aggregating ", resolution, " exposure: ", basename(file))
    arrow::open_dataset(file) |>
      dplyr::select(dplyr::all_of(c(groups, value_col))) |>
      dplyr::rename(value = dplyr::all_of(value_col)) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(groups))) |>
      dplyr::summarise(n_days = sum(!is.na(value)),
                       value = mean(value, na.rm = TRUE), .groups = "drop") |>
      dplyr::filter(n_days >= minimum_days, is.finite(value)) |>
      dplyr::collect()
  })
  result <- dplyr::bind_rows(pieces) |>
    dplyr::select(-n_days) |>
    dplyr::rename(!!out_col := value) |>
    dplyr::mutate(zip = sprintf("%05d", as.integer(zip)), year = as.integer(year))
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(result, cache_file)
  result
}
