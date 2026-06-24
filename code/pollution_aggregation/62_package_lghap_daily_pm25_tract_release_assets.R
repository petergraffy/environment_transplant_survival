#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- file.path("..", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- file.path("..", "..", "code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}
local_appdata <- Sys.getenv("LOCALAPPDATA", unset = NA_character_)
if (!is.na(local_appdata) && nzchar(local_appdata)) {
  r_minor <- strsplit(R.version$minor, "[.]", fixed = FALSE)[[1]][[1]]
  user_lib <- file.path(local_appdata, "R", "win-library", paste0(R.version$major, ".", r_minor))
  if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
}

region_slug <- Sys.getenv("TRACT_PM25_REGION_SLUG", "cook_county_il")
region_label <- Sys.getenv("TRACT_PM25_REGION_LABEL", "Cook County, Illinois")
source_dir <- Sys.getenv("TRACT_PM25_SOURCE_DIR", file.path("data", "processed", paste0("lghap_pm25_tract_daily_", region_slug)))
release_dir <- Sys.getenv("TRACT_PM25_RELEASE_DIR", file.path("data", "release", paste0("lghap_pm25_tract_daily_", region_slug, "_csv")))
dir.create(release_dir, recursive = TRUE, showWarnings = FALSE)

parse_int_set <- function(value, default) {
  if (!nzchar(value)) return(default)
  pieces <- strsplit(value, ",")[[1]]
  out <- integer()
  for (piece in trimws(pieces)) {
    bounds <- as.integer(strsplit(piece, ":", fixed = TRUE)[[1]])
    if (length(bounds) == 2L) {
      out <- c(out, seq(bounds[1], bounds[2]))
    } else {
      out <- c(out, bounds[1])
    }
  }
  unique(out)
}

years <- parse_int_set(Sys.getenv("TRACT_PM25_RELEASE_YEARS"), 2005:2021)
overwrite <- identical(tolower(Sys.getenv("TRACT_PM25_RELEASE_OVERWRITE")), "true")

value_cols <- c(
  "tract_geoid", "date", "pollutant", "year", "month", "temporal_resolution",
  "pm25_ug_m3", "value_source", "fill_distance_m", "fill_cell"
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
  files <- file.path(source_dir, sprintf("lghap_pm25_tract_daily_%04d_%02d.csv.gz", year, 1:12))
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Missing source files: ", paste(missing, collapse = ", "), call. = FALSE)

  do.call(rbind, lapply(files, function(path) {
    log_msg("Reading ", basename(path))
    dat <- read.csv(
      gzfile(path),
      colClasses = c(
        tract_geoid = "character",
        date = "character",
        pollutant = "character",
        year = "integer",
        month = "integer",
        temporal_resolution = "character",
        pm25_ug_m3 = "numeric",
        value_source = "character",
        fill_distance_m = "numeric",
        fill_cell = "numeric"
      )
    )
    dat[, value_cols]
  }))
}

convert_year <- function(year) {
  dest <- file.path(release_dir, sprintf("lghap_pm25_tract_daily_%04d.csv.gz", year))

  if (!file.exists(dest) || overwrite) {
    dat <- read_year(year)
    dat$tract_geoid <- sprintf("%011s", dat$tract_geoid)
    dat$date <- as.Date(dat$date)
    dat$year <- as.integer(dat$year)
    dat$month <- as.integer(dat$month)
    log_msg("Writing CSV.gz for ", year)
    con <- gzfile(dest, open = "wt")
    on.exit(close(con), add = TRUE)
    write.csv(dat, con, row.names = FALSE, na = "")
  } else {
    log_msg("Using existing CSV.gz for ", year)
    dat <- read_year(year)
  }

  info <- file.info(dest)
  data.frame(
    year = year,
    file = basename(dest),
    path = normalizePath(dest, winslash = "/", mustWork = TRUE),
    rows = nrow(dat),
    expected_rows = length(unique(dat$tract_geoid)) * ifelse(is_leap_year(year), 366L, 365L),
    tracts = length(unique(dat$tract_geoid)),
    dates = length(unique(dat$date)),
    first_date = min(dat$date, na.rm = TRUE),
    last_date = max(dat$date, na.rm = TRUE),
    missing_values = sum(is.na(dat$pm25_ug_m3)),
    mean_pm25_ug_m3 = mean(dat$pm25_ug_m3, na.rm = TRUE),
    size_bytes = as.numeric(info$size),
    size_mb = round(as.numeric(info$size) / 1024^2, 1),
    sha256 = sha256_file(dest)
  )
}

log_msg("Packaging LGHAP daily PM2.5 tract CSV.gz release assets")
manifest <- do.call(rbind, lapply(years, convert_year))
manifest <- manifest[order(manifest$year), ]
write.csv(manifest, file.path(release_dir, "lghap_pm25_tract_daily_csv_manifest.csv"), row.names = FALSE)

qc_all_path <- file.path(source_dir, "lghap_pm25_tract_daily_qc_all_months.csv")
if (file.exists(qc_all_path)) {
  qc_all <- read.csv(qc_all_path)
  qc_all <- qc_all[qc_all$year %in% years, ]
  write.csv(qc_all, file.path(release_dir, "lghap_pm25_tract_daily_qc_all_months.csv"), row.names = FALSE)
}

source_manifest <- file.path(source_dir, "lghap_pm25_tract_daily_manifest.csv")
if (file.exists(source_manifest)) {
  src <- read.csv(source_manifest)
  src <- src[src$year %in% years, ]
  write.csv(src, file.path(release_dir, "lghap_pm25_tract_daily_source_manifest.csv"), row.names = FALSE)
}

readme <- c(
  paste0("# LGHAP Daily PM2.5 Census Tract CSV.gz Release Assets: ", region_label),
  "",
  paste0("These files contain daily LGHAP PM2.5 exposures aggregated to 2020 census tract polygons for ", region_label, "."),
  "",
  "Assets:",
  "",
  "- `lghap_pm25_tract_daily_YYYY.csv.gz`: yearly daily tract PM2.5 values, one row per tract-date.",
  "- `lghap_pm25_tract_daily_csv_manifest.csv`: row counts, date coverage, completeness, file sizes, and SHA-256 checksums.",
  "- `lghap_pm25_tract_daily_qc_all_months.csv`: monthly completeness and distribution QC from the source aggregation.",
  "- `lghap_pm25_tract_daily_source_manifest.csv`: monthly source CSV manifest, when available.",
  "- `README.md`: variables and citation notes.",
  "",
  "Columns:",
  "",
  paste0("- `", value_cols, "`"),
  "",
  "The `tract_geoid` column is the 11-digit 2020 census tract GEOID. The `pm25_ug_m3` column is daily PM2.5 in micrograms per cubic meter. `value_source`, `fill_distance_m`, and `fill_cell` provide audit metadata for nearest-raster-cell fills used when polygon zonal means are unavailable.",
  "",
  "Citation:",
  "",
  "Bai, K., Li, K., Shao, L., Li, X., Liu, C., Li, Z., Ma, M., Han, D., Sun, Y., Zheng, Z., Li, R., Chang, N.-B., and Guo, J.: LGHAP v2: a global gap-free aerosol optical depth and PM2.5 concentration dataset since 2000 derived via big Earth data analytics, Earth Syst. Sci. Data, 16, 2425-2448, https://doi.org/10.5194/essd-16-2425-2024, 2024.",
  "",
  "Also cite this repository/release for the census tract aggregation workflow. The reproducible aggregation script is `code/61_aggregate_lghap_daily_pm25_to_tract.R` and this packaging script is `code/pollution_aggregation/62_package_lghap_daily_pm25_tract_release_assets.R`."
)
writeLines(readme, file.path(release_dir, "README.md"), useBytes = TRUE)

log_msg("Wrote release assets to ", normalizePath(release_dir, winslash = "/"))
print(data.frame(
  years = nrow(manifest),
  total_size_mb = sum(manifest$size_mb),
  max_asset_mb = max(manifest$size_mb),
  missing_values = sum(manifest$missing_values)
))
