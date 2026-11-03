#!/usr/bin/env Rscript

region_slug <- Sys.getenv("TRACT_O3_REGION_SLUG", "cook_county_il")
region_label <- Sys.getenv("TRACT_O3_REGION_LABEL", "Cook County, Illinois")
source_dir <- Sys.getenv("TRACT_O3_SOURCE_DIR", file.path("data", "processed", paste0("o3_tract_daily_", region_slug)))
release_dir <- Sys.getenv("TRACT_O3_RELEASE_DIR", file.path("data", "release", paste0("o3_tract_daily_", region_slug, "_csv")))
dir.create(release_dir, recursive = TRUE, showWarnings = FALSE)

parse_int_set <- function(value, default) {
  if (!nzchar(value)) return(default)
  pieces <- strsplit(value, ",")[[1]]
  out <- integer()
  for (piece in trimws(pieces)) {
    bounds <- as.integer(strsplit(piece, ":", fixed = TRUE)[[1]])
    if (length(bounds) == 2L) out <- c(out, seq(bounds[1], bounds[2])) else out <- c(out, bounds[1])
  }
  unique(out)
}

years <- parse_int_set(Sys.getenv("TRACT_O3_RELEASE_YEARS"), 2005:2024)
overwrite <- identical(tolower(Sys.getenv("TRACT_O3_RELEASE_OVERWRITE")), "true")
value_cols <- c("tract_geoid", "date", "pollutant", "year", "month", "temporal_resolution", "o3_ppb", "value_source", "fill_distance_m", "fill_cell")

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

sha256_file <- function(path) unname(tools::sha256sum(path))
is_leap_year <- function(year) (year %% 4L == 0L & year %% 100L != 0L) | year %% 400L == 0L

read_year <- function(year) {
  files <- file.path(source_dir, sprintf("o3_tract_daily_%04d_%02d.csv.gz", year, 1:12))
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
        o3_ppb = "numeric",
        value_source = "character",
        fill_distance_m = "numeric",
        fill_cell = "numeric"
      )
    )
    dat[, value_cols]
  }))
}

convert_year <- function(year) {
  dest <- file.path(release_dir, sprintf("o3_tract_daily_%04d.csv.gz", year))
  if (!file.exists(dest) || overwrite) {
    dat <- read_year(year)
    dat$tract_geoid <- sprintf("%011s", dat$tract_geoid)
    dat$date <- as.Date(dat$date)
    dat$year <- as.integer(dat$year)
    dat$month <- as.integer(dat$month)
    log_msg("Writing CSV.gz for ", year)
    con <- gzfile(dest, open = "wt")
    write.csv(dat, con, row.names = FALSE, na = "")
    close(con)
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
    missing_values = sum(is.na(dat$o3_ppb)),
    mean_o3_ppb = mean(dat$o3_ppb, na.rm = TRUE),
    size_bytes = as.numeric(info$size),
    size_mb = round(as.numeric(info$size) / 1024^2, 1),
    sha256 = sha256_file(dest)
  )
}

log_msg("Packaging daily ozone tract CSV.gz release assets")
manifest <- do.call(rbind, lapply(years, convert_year))
manifest <- manifest[order(manifest$year), ]
write.csv(manifest, file.path(release_dir, "o3_tract_daily_csv_manifest.csv"), row.names = FALSE)

qc_all_path <- file.path(source_dir, "o3_tract_daily_qc_all_months.csv")
if (file.exists(qc_all_path)) {
  qc_all <- read.csv(qc_all_path)
  qc_all <- qc_all[qc_all$year %in% years, ]
  write.csv(qc_all, file.path(release_dir, "o3_tract_daily_qc_all_months.csv"), row.names = FALSE)
}

source_manifest <- file.path(source_dir, "o3_tract_daily_manifest.csv")
if (file.exists(source_manifest)) {
  src <- read.csv(source_manifest)
  src <- src[src$year %in% years, ]
  write.csv(src, file.path(release_dir, "o3_tract_daily_source_manifest.csv"), row.names = FALSE)
}

readme <- c(
  paste0("# Daily Ozone Census Tract CSV.gz Release Assets: ", region_label),
  "",
  paste0("These files contain daily MDA8 ozone exposures aggregated to 2020 census tract polygons for ", region_label, "."),
  "",
  "Assets:",
  "",
  "- `o3_tract_daily_YYYY.csv.gz`: yearly daily tract ozone values, one row per tract-date.",
  "- `o3_tract_daily_csv_manifest.csv`: row counts, date coverage, completeness, file sizes, and SHA-256 checksums.",
  "- `o3_tract_daily_qc_all_months.csv`: monthly completeness and distribution QC from the source aggregation.",
  "- `o3_tract_daily_source_manifest.csv`: monthly source CSV manifest, when available.",
  "- `README_o3.md`: variables and citation notes.",
  "",
  "Columns:",
  "",
  paste0("- `", value_cols, "`"),
  "",
  "The `tract_geoid` column is the 11-digit 2020 census tract GEOID. The `o3_ppb` column is daily maximum 8-hour average ozone in parts per billion. `value_source`, `fill_distance_m`, and `fill_cell` provide audit metadata for nearest-raster-cell fills used when polygon zonal means are unavailable.",
  "",
  "Citation:",
  "",
  "Please cite the source ozone estimates using the Science article DOI supplied with these data: https://www.science.org/doi/10.1126/science.aed3197.",
  "",
  "Also cite this repository/release for the census tract aggregation workflow. The reproducible aggregation script is `code/64_aggregate_daily_o3_to_tract.R` and this packaging script is `code/pollution_aggregation/65_package_daily_o3_tract_release_assets.R`."
)
writeLines(readme, file.path(release_dir, "README_o3.md"), useBytes = TRUE)

log_msg("Wrote release assets to ", normalizePath(release_dir, winslash = "/"))
print(data.frame(
  years = nrow(manifest),
  total_size_mb = sum(manifest$size_mb),
  max_asset_mb = max(manifest$size_mb),
  missing_values = sum(manifest$missing_values)
))
