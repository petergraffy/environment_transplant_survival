#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(readr)
})

source_dir <- file.path("output", "zip_pollution", "filled")
release_dir <- file.path("data", "release", "air_pollution_zcta_parquet")
dir.create(release_dir, recursive = TRUE, showWarnings = FALSE)

assets <- tibble::tribble(
  ~pollutant, ~source_file, ~value_col, ~output_file,
  "pm25", "all_pm25_zip_filled.csv.gz", "pm25_ug_m3", "air_pollution_zcta_pm25_monthly_2005_2023.parquet",
  "o3", "all_o3_zip_filled.csv.gz", "o3_ppb", "air_pollution_zcta_o3_monthly_2005_2023.parquet",
  "no2", "all_no2_zip_filled.csv.gz", "no2", "air_pollution_zcta_no2_annual_2005_2024.parquet"
)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

sha256_file <- function(path) {
  unname(tools::sha256sum(path))
}

convert_asset <- function(pollutant, source_file, value_col, output_file) {
  src <- file.path(source_dir, source_file)
  dest <- file.path(release_dir, output_file)
  if (!file.exists(src)) stop("Missing source file: ", src, call. = FALSE)

  if (!file.exists(dest) || identical(tolower(Sys.getenv("AIR_POLLUTION_RELEASE_OVERWRITE")), "true")) {
    log_msg("Reading ", source_file)
    dat <- read_csv(src, col_types = cols(zip = col_character(), .default = col_guess()), show_col_types = FALSE)
    dat$zip <- sprintf("%05s", dat$zip)
    dat$year <- as.integer(dat$year)
    dat$month <- as.integer(dat$month)
    log_msg("Writing ", output_file)
    write_parquet(
      dat,
      dest,
      compression = "zstd",
      compression_level = 9,
      use_dictionary = TRUE
    )
  } else {
    log_msg("Using existing ", output_file)
    dat <- read_parquet(dest, as_data_frame = TRUE)
  }

  info <- file.info(dest)
  tibble::tibble(
    pollutant = pollutant,
    file = basename(dest),
    path = normalizePath(dest, winslash = "/", mustWork = TRUE),
    rows = nrow(dat),
    zips = dplyr::n_distinct(dat$zip),
    first_year = min(dat$year, na.rm = TRUE),
    last_year = max(dat$year, na.rm = TRUE),
    temporal_resolution = paste(sort(unique(dat$temporal_resolution)), collapse = ";"),
    value_column = value_col,
    missing_values = sum(is.na(dat[[value_col]])),
    size_bytes = as.numeric(info$size),
    size_mb = round(as.numeric(info$size) / 1024^2, 1),
    sha256 = sha256_file(dest)
  )
}

manifest <- bind_rows(Map(
  convert_asset,
  assets$pollutant,
  assets$source_file,
  assets$value_col,
  assets$output_file
))

write_csv(manifest, file.path(release_dir, "air_pollution_zcta_parquet_manifest.csv"))

meta_files <- c(
  "fill_summary.csv",
  file.path("..", "fill_audit_summary.csv")
)
for (rel in meta_files) {
  src <- normalizePath(file.path(source_dir, rel), winslash = "/", mustWork = FALSE)
  if (file.exists(src)) file.copy(src, file.path(release_dir, basename(src)), overwrite = TRUE)
}

readme <- c(
  "# Air Pollution ZCTA Parquet Release Assets",
  "",
  "These files contain final filled air-pollution exposure surfaces aggregated to 2020 ZCTA5 polygons.",
  "",
  "Assets:",
  "",
  "- `air_pollution_zcta_pm25_monthly_2005_2023.parquet`: monthly PM2.5, `pm25_ug_m3`.",
  "- `air_pollution_zcta_o3_monthly_2005_2023.parquet`: monthly ozone, `o3_ppb`.",
  "- `air_pollution_zcta_no2_annual_2005_2024.parquet`: annual NO2, `no2`.",
  "- `air_pollution_zcta_parquet_manifest.csv`: sizes, row counts, completeness, and SHA-256 checksums.",
  "- `fill_summary.csv`: nearest-fill audit by source file.",
  "- `fill_audit_summary.csv`: additional fill audit summary when available.",
  "",
  "Columns include `zip`, `pollutant`, `year`, `month`, `temporal_resolution`, the pollutant value column, `value_source`, and `fill_distance_km`.",
  "PM2.5 and ozone currently cover 2005-2023 in the local source files; NO2 covers 2005-2024."
)
writeLines(readme, file.path(release_dir, "README.md"), useBytes = TRUE)

log_msg("Wrote release assets to ", normalizePath(release_dir, winslash = "/"))
print(manifest %>% summarise(files = n(), total_size_mb = sum(size_mb), max_asset_mb = max(size_mb), missing_values = sum(missing_values)))
