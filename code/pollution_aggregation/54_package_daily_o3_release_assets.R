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

source_dir <- file.path("data", "processed", "o3_zcta_daily")
release_dir <- file.path("data", "release", "o3_zcta_daily_parquet")
dir.create(release_dir, recursive = TRUE, showWarnings = FALSE)

years <- as.integer(strsplit(Sys.getenv("DAILY_O3_RELEASE_YEARS", "2005:2024"), ":", fixed = TRUE)[[1]])
if (length(years) == 2L) years <- seq(years[1], years[2])

value_cols <- c(
  "zip", "date", "pollutant", "year", "month", "temporal_resolution",
  "o3_ppb", "value_source", "fill_distance_m", "fill_cell"
)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

sha256_file <- function(path) {
  unname(tools::sha256sum(path))
}

is_leap_year <- function(year) {
  (year %% 4L == 0L & year %% 100L != 0L) | year %% 400L == 0L
}

read_year <- function(year) {
  files <- file.path(source_dir, sprintf("o3_zcta_daily_%04d_%02d.csv.gz", year, 1:12))
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Missing source files: ", paste(missing, collapse = ", "), call. = FALSE)

  bind_rows(lapply(files, function(path) {
    log_msg("Reading ", basename(path))
    read_csv(
      path,
      col_select = all_of(value_cols),
      col_types = cols(
        zip = col_character(),
        date = col_date(),
        pollutant = col_character(),
        year = col_integer(),
        month = col_integer(),
        temporal_resolution = col_character(),
        o3_ppb = col_double(),
        value_source = col_character(),
        fill_distance_m = col_double(),
        fill_cell = col_double()
      ),
      show_col_types = FALSE,
      progress = FALSE
    )
  }))
}

convert_year <- function(year) {
  dest <- file.path(release_dir, sprintf("o3_zcta_daily_%04d.parquet", year))

  if (!file.exists(dest) || identical(tolower(Sys.getenv("DAILY_O3_RELEASE_OVERWRITE")), "true")) {
    dat <- read_year(year)
    dat$zip <- sprintf("%05s", dat$zip)
    dat$year <- as.integer(dat$year)
    dat$month <- as.integer(dat$month)
    log_msg("Writing Parquet for ", year)
    write_parquet(
      dat,
      dest,
      compression = "zstd",
      compression_level = 9,
      use_dictionary = TRUE
    )
  } else {
    log_msg("Using existing Parquet for ", year)
    dat <- read_parquet(dest, as_data_frame = TRUE)
  }

  info <- file.info(dest)
  tibble::tibble(
    year = year,
    file = basename(dest),
    path = normalizePath(dest, winslash = "/", mustWork = TRUE),
    rows = nrow(dat),
    expected_rows = 33300L * ifelse(is_leap_year(year), 366L, 365L),
    zips = dplyr::n_distinct(dat$zip),
    dates = dplyr::n_distinct(dat$date),
    first_date = min(dat$date, na.rm = TRUE),
    last_date = max(dat$date, na.rm = TRUE),
    missing_values = sum(is.na(dat$o3_ppb)),
    mean_o3_ppb = mean(dat$o3_ppb, na.rm = TRUE),
    size_bytes = as.numeric(info$size),
    size_mb = round(as.numeric(info$size) / 1024^2, 1),
    sha256 = sha256_file(dest)
  )
}

log_msg("Packaging daily ozone ZCTA Parquet release assets")
manifest <- bind_rows(lapply(years, convert_year))
write_csv(manifest, file.path(release_dir, "o3_zcta_daily_parquet_manifest.csv"))

qc_files <- file.path(source_dir, sprintf("o3_zcta_daily_qc_%04d_all_months.csv", years))
qc_files <- qc_files[file.exists(qc_files)]
if (length(qc_files)) {
  qc_all <- bind_rows(lapply(qc_files, read_csv, show_col_types = FALSE))
  write_csv(qc_all, file.path(release_dir, "o3_zcta_daily_qc_all_years.csv"))
}

readme <- c(
  "# Daily Ozone ZCTA Parquet Release Assets",
  "",
  "These files contain daily MDA8 ozone exposures aggregated to 2020 ZCTA5 polygons for CONUS.",
  "",
  "Assets:",
  "",
  "- `o3_zcta_daily_YYYY.parquet`: yearly daily ZCTA ozone values, one row per ZCTA-date.",
  "- `o3_zcta_daily_parquet_manifest.csv`: row counts, date coverage, completeness, file sizes, and SHA-256 checksums.",
  "- `o3_zcta_daily_qc_all_years.csv`: monthly completeness and distribution QC from the source aggregation.",
  "- `README.md`: variables and citation notes.",
  "",
  "Columns:",
  "",
  paste0("- `", value_cols, "`"),
  "",
  "The `o3_ppb` column is daily maximum 8-hour average ozone in parts per billion. `value_source`, `fill_distance_m`, and `fill_cell` provide audit metadata for nearest-cell fills; the completed 2005-2024 local outputs have zero missing ZCTA-day ozone values.",
  "",
  "Citation:",
  "",
  "Please cite the source ozone estimates using the Science article DOI supplied with these data: https://www.science.org/doi/10.1126/science.aed3197.",
  "",
  "Also cite this repository/release for the ZCTA aggregation workflow. The reproducible aggregation script is `code/pollution_aggregation/52_aggregate_daily_o3_to_zcta.R`, the QA mapping script is `code/pollution_aggregation/53_map_daily_o3_zcta_examples.R`, and this packaging script is `code/pollution_aggregation/54_package_daily_o3_release_assets.R`."
)
writeLines(readme, file.path(release_dir, "README.md"), useBytes = TRUE)

log_msg("Wrote release assets to ", normalizePath(release_dir, winslash = "/"))
print(manifest %>% summarise(
  years = n(),
  total_size_mb = sum(size_mb),
  max_asset_mb = max(size_mb),
  missing_values = sum(missing_values)
))
