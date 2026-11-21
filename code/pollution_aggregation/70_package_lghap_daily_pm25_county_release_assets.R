#!/usr/bin/env Rscript
source(file.path("code", "r_runtime.R"))
ensure_user_library()
suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

source_dir <- "data/processed/lghap_pm25_county_daily"
release_dir <- "data/release/lghap_pm25_county_parquet"
dir.create(release_dir, recursive = TRUE, showWarnings = FALSE)
years <- 2005:2024
fips_types <- c(county_fips = "character", state_fips = "character")
qc <- fread(file.path(source_dir, "lghap_pm25_county_daily_manifest.csv"))
annual <- fread(file.path(source_dir, "lghap_pm25_county_annual.csv"), colClasses = fips_types)
lookup <- fread(file.path(source_dir, "county_lookup_2020_conus.csv"), colClasses = fips_types)
fill_audit <- fread(file.path(source_dir, "lghap_pm25_county_fill_audit.csv"), colClasses = fips_types)
stopifnot(nrow(qc) == 240L, !anyDuplicated(qc[, .(year, month)]), setequal(qc$year, years),
  all(qc$month %in% 1:12), all(qc$missing_values == 0L),
  nrow(lookup) == 3108L, !anyDuplicated(lookup$county_fips),
  all(grepl("^[0-9]{5}$", lookup$county_fips)),
  nrow(annual) == nrow(lookup) * length(years), setequal(annual$year, years),
  !anyDuplicated(annual[, .(county_fips, year)]),
  nrow(fill_audit) == sum(qc$nearest_fill_values),
  !anyDuplicated(fill_audit[, .(county_fips, date)]),
  identical(unname(tools::sha256sum(file.path(source_dir, qc$file))), qc$sha256))

read_month <- function(file) {
  con <- gzfile(file.path(source_dir, file), "rt")
  on.exit(close(con))
  fread(text = paste(readLines(con), collapse = "\n"), colClasses = fips_types)
}

write_asset <- function(dat, file, resolution, first_date, last_date, fills, year = NA_integer_) {
  path <- file.path(release_dir, file)
  stage <- tempfile(tmpdir = release_dir)
  on.exit(unlink(stage))
  write_parquet(dat, stage, compression = "zstd", compression_level = 9, use_dictionary = TRUE)
  reopened <- as.data.table(read_parquet(stage))
  stopifnot(isTRUE(all.equal(reopened, dat, tolerance = 0, check.attributes = FALSE)),
    is.character(reopened$county_fips), is.character(reopened$state_fips),
    is.numeric(reopened$pm25_ug_m3),
    resolution != "daily" || inherits(reopened$date, "Date"))
  if (!file.rename(stage, path)) stop("Could not finalize ", path)
  data.table(file = file, temporal_resolution = resolution, year = year, rows = nrow(dat),
    counties = uniqueN(dat$county_fips), first_date = first_date, last_date = last_date,
    missing_values = sum(!is.finite(dat$pm25_ug_m3)), nearest_fill_days = fills,
    size_bytes = file.size(path), sha256 = unname(tools::sha256sum(path)))
}

daily_manifest <- rbindlist(lapply(years, function(y) {
  message("Packaging and checking ", y)
  months <- qc[year == y][order(month)]
  stopifnot(identical(months$month, 1:12))
  dat <- rbindlist(lapply(months$file, read_month))
  dat[, `:=`(date = as.Date(date), pm25_ug_m3 = as.double(pm25_ug_m3),
    valid_area_fraction = as.double(valid_area_fraction),
    fill_distance_m = as.double(fill_distance_m), fill_cell = as.integer(fill_cell))]
  dates <- seq(as.Date(sprintf("%d-01-01", y)), as.Date(sprintf("%d-12-31", y)), by = "day")
  stopifnot(setequal(dat$county_fips, lookup$county_fips), setequal(dat$date, dates),
    all(dat$year == y), nrow(dat) == length(dates) * nrow(lookup),
    !anyDuplicated(dat[, .(county_fips, date)]), all(is.finite(dat$pm25_ug_m3)),
    all(dat$pm25_ug_m3 >= 0), all(dat$state_fips == substr(dat$county_fips, 1L, 2L)))

  daily_fills <- dat[value_source == "nearest_raster_cell", .(county_fips, date, fill_distance_m, fill_cell)]
  audit_fills <- fill_audit[year == y, .(county_fips, date = as.Date(date), fill_distance_m, fill_cell)]
  stopifnot(nrow(daily_fills) == sum(months$nearest_fill_values),
    isTRUE(all.equal(daily_fills, audit_fills, check.attributes = FALSE)))

  calculated <- dat[, .(days = .N, pm25_ug_m3 = mean(pm25_ug_m3),
    nearest_fill_days = sum(value_source == "nearest_raster_cell")), by = county_fips]
  expected <- annual[year == y, .(county_fips, days, pm25_ug_m3, nearest_fill_days)]
  stopifnot(isTRUE(all.equal(calculated[order(county_fips)], expected[order(county_fips)],
    tolerance = 1e-12, check.attributes = FALSE)))
  write_asset(dat, sprintf("lghap_pm25_county_daily_%d.parquet", y), "daily",
    min(dates), max(dates), nrow(daily_fills), y)
}))

annual_manifest <- write_asset(annual, "lghap_pm25_county_annual.parquet", "annual",
  min(daily_manifest$first_date), max(daily_manifest$last_date), sum(annual$nearest_fill_days))
manifest <- rbind(daily_manifest, annual_manifest)
fwrite(manifest, file.path(release_dir, "lghap_pm25_county_parquet_manifest.csv"))
metadata <- c("lghap_pm25_county_daily_manifest.csv", "lghap_pm25_county_fill_audit.csv",
  "county_lookup_2020_conus.csv")
stopifnot(all(file.copy(file.path(source_dir, metadata), release_dir, overwrite = TRUE)),
  file.copy("docs/lghap_daily_pm25_county.md", file.path(release_dir, "README.md"), overwrite = TRUE))
message("Verified ", nrow(manifest), " Parquet assets in ", release_dir)
