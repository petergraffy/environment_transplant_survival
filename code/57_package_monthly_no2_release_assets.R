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

source_dir <- file.path("data", "processed", "no2_zcta_monthly")
release_dir <- file.path("data", "release", "no2_zcta_monthly_parquet")
dir.create(release_dir, recursive = TRUE, showWarnings = FALSE)

years <- as.integer(strsplit(Sys.getenv("MONTHLY_NO2_RELEASE_YEARS", "2019:2025"), ":", fixed = TRUE)[[1]])
if (length(years) == 2L) years <- seq(years[1], years[2])

value_cols <- c(
  "zip", "pollutant", "year", "month", "temporal_resolution",
  "no2_ppbv", "value_source", "fill_distance_m", "fill_cell"
)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

sha256_file <- function(path) {
  unname(tools::sha256sum(path))
}

read_year <- function(year) {
  files <- file.path(source_dir, sprintf("no2_zcta_monthly_%04d_%02d.csv.gz", year, 1:12))
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Missing source files: ", paste(missing, collapse = ", "), call. = FALSE)

  bind_rows(lapply(files, function(path) {
    log_msg("Reading ", basename(path))
    read_csv(
      path,
      col_select = all_of(value_cols),
      col_types = cols(
        zip = col_character(),
        pollutant = col_character(),
        year = col_integer(),
        month = col_integer(),
        temporal_resolution = col_character(),
        no2_ppbv = col_double(),
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
  dest <- file.path(release_dir, sprintf("no2_zcta_monthly_%04d.parquet", year))

  if (!file.exists(dest) || identical(tolower(Sys.getenv("MONTHLY_NO2_RELEASE_OVERWRITE")), "true")) {
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
    expected_rows = 33300L * 12L,
    zips = dplyr::n_distinct(dat$zip),
    months = dplyr::n_distinct(dat$month),
    first_month = min(dat$month, na.rm = TRUE),
    last_month = max(dat$month, na.rm = TRUE),
    missing_values = sum(is.na(dat$no2_ppbv)),
    mean_no2_ppbv = mean(dat$no2_ppbv, na.rm = TRUE),
    size_bytes = as.numeric(info$size),
    size_mb = round(as.numeric(info$size) / 1024^2, 1),
    sha256 = sha256_file(dest)
  )
}

log_msg("Packaging monthly NO2 ZCTA Parquet release assets")
manifest <- bind_rows(lapply(years, convert_year))
write_csv(manifest, file.path(release_dir, "no2_zcta_monthly_parquet_manifest.csv"))

qc_src <- file.path(source_dir, "no2_zcta_monthly_qc_all_months.csv")
if (file.exists(qc_src)) {
  file.copy(qc_src, file.path(release_dir, basename(qc_src)), overwrite = TRUE)
}

readme <- c(
  "# Monthly NO2 ZCTA Parquet Release Assets",
  "",
  "These files contain monthly TROPOMI/LUR CONUS surface NO2 estimates aggregated to 2020 ZCTA5 polygons.",
  "",
  "Assets:",
  "",
  "- `no2_zcta_monthly_YYYY.parquet`: yearly monthly ZCTA NO2 values, one row per ZCTA-month.",
  "- `no2_zcta_monthly_parquet_manifest.csv`: row counts, month coverage, completeness, file sizes, and SHA-256 checksums.",
  "- `no2_zcta_monthly_qc_all_months.csv`: monthly completeness and distribution QC from the source aggregation.",
  "- `README.md`: variables and source notes.",
  "",
  "Columns:",
  "",
  paste0("- `", value_cols, "`"),
  "",
  "The `no2_ppbv` column is monthly surface NO2 in parts per billion by volume. `value_source`, `fill_distance_m`, and `fill_cell` provide audit metadata for nearest-cell fills; the completed 2019-2025 local outputs have zero missing ZCTA-month NO2 values.",
  "",
  "Citation:",
  "",
  "Nawaz, M. Omar. Monthly and Annual US TROPOMI Surface NO2 Estimates (~1km x 1km), Version v2. Zenodo, 2026. https://doi.org/10.5281/zenodo.18919769.",
  "",
  "Zenodo record: https://zenodo.org/records/18919769.",
  "Related publication DOI listed by Zenodo: https://doi.org/10.1021/acsestair.4c00153.",
  "Also cite this repository/release for the ZCTA aggregation workflow.",
  "",
  "The reproducible aggregation script is `code/56_aggregate_monthly_no2_to_zcta.R` and this packaging script is `code/57_package_monthly_no2_release_assets.R`."
)
writeLines(readme, file.path(release_dir, "README.md"), useBytes = TRUE)

log_msg("Wrote release assets to ", normalizePath(release_dir, winslash = "/"))
print(manifest %>% summarise(
  years = n(),
  total_size_mb = sum(size_mb),
  max_asset_mb = max(size_mb),
  missing_values = sum(missing_values)
))
