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

source_dir <- file.path("data", "processed", "gridmet_zcta_daily")
release_dir <- file.path("data", "release", "gridmet_zcta_daily_parquet")
dir.create(release_dir, recursive = TRUE, showWarnings = FALSE)

years <- as.integer(strsplit(Sys.getenv("GRIDMET_RELEASE_YEARS", "2005:2025"), ":", fixed = TRUE)[[1]])
if (length(years) == 2L) years <- seq(years[1], years[2])

value_cols <- c(
  "zip", "date", "year",
  "tmax_c", "tmin_c", "tmean_c", "temp_range_c",
  "rmax_pct", "rmin_pct", "rhmean_pct",
  "pr", "sph", "srad", "vs", "pet", "etr", "erc"
)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

sha256_file <- function(path) {
  unname(tools::sha256sum(path))
}

convert_year <- function(year) {
  src <- file.path(source_dir, sprintf("gridmet_zcta_daily_%04d.csv.gz", year))
  dest <- file.path(release_dir, sprintf("gridmet_zcta_daily_%04d.parquet", year))
  if (!file.exists(src)) stop("Missing source file: ", src, call. = FALSE)

  if (!file.exists(dest) || identical(tolower(Sys.getenv("GRIDMET_RELEASE_OVERWRITE")), "true")) {
    log_msg("Reading value columns for ", year)
    dat <- read_csv_arrow(src, col_select = all_of(value_cols), as_data_frame = TRUE)
    dat$zip <- sprintf("%05d", as.integer(dat$zip))
    dat$year <- as.integer(dat$year)
    log_msg("Writing Parquet for ", year)
    write_parquet(
      dat,
      dest,
      compression = "zstd",
      compression_level = 9,
      use_dictionary = TRUE
    )
    rm(dat)
    gc()
  } else {
    log_msg("Using existing Parquet for ", year)
  }

  info <- file.info(dest)
  tibble::tibble(
    year = year,
    file = basename(dest),
    path = normalizePath(dest, winslash = "/", mustWork = TRUE),
    rows = ifelse(year %% 4L == 0L, 33300L * 366L, 33300L * 365L),
    size_bytes = as.numeric(info$size),
    size_mb = round(as.numeric(info$size) / 1024^2, 1),
    sha256 = sha256_file(dest)
  )
}

log_msg("Packaging gridMET ZCTA daily Parquet release assets")
manifest <- bind_rows(lapply(years, convert_year))
write_csv(manifest, file.path(release_dir, "gridmet_zcta_daily_parquet_manifest.csv"))

qc_src <- file.path(source_dir, "gridmet_zcta_daily_qc_all_years.csv")
if (file.exists(qc_src)) {
  file.copy(qc_src, file.path(release_dir, basename(qc_src)), overwrite = TRUE)
}

readme <- c(
  "# gridMET ZCTA Daily Parquet Release Assets",
  "",
  "These files contain value-only daily gridMET exposures aggregated to 2020 ZCTA5 polygons.",
  "",
  "Each yearly Parquet file contains:",
  "",
  paste0("- `", value_cols, "`"),
  "",
  "Audit/fill summaries are in `gridmet_zcta_daily_qc_all_years.csv` and the committed documentation.",
  "The original local processing outputs include source/fill columns and can be regenerated with `code/exploratory/20_aggregate_gridmet_to_zcta.R`."
)
writeLines(readme, file.path(release_dir, "README.md"), useBytes = TRUE)

log_msg("Wrote release assets to ", normalizePath(release_dir, winslash = "/"))
print(manifest %>% summarise(years = n(), total_size_mb = sum(size_mb), max_asset_mb = max(size_mb)))
