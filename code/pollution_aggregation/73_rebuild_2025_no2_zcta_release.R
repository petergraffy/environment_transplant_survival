#!/usr/bin/env Rscript
source("code/r_runtime.R")
ensure_user_library()
source("code/pollution_aggregation/aggregate_pollution_to_zcta.R")
source("code/pollution_aggregation/fill_missing_zip_pollution_values.R")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, file.exists(args[1]), grepl("_2025\\.v.*\\.nc$", basename(args[1])))
input <- args[1]
release_dir <- "data/release/air_pollution_zcta_parquet"
asset <- file.path(release_dir, "air_pollution_zcta_no2_annual_2005_2025.parquet")
period <- "no2_2025.csv.gz"

raw <- aggregate_no2_period(input, read_zctas(), 2025L)
filled <- fill_one(raw)
old <- as.data.table(arrow::read_parquet(asset))
historical <- old[year != 2025L]
stopifnot(nrow(filled) == 33300L, all(grepl("^[0-9]{5}$", filled$zip)),
  !anyDuplicated(filled$zip), setequal(filled$zip, old[year == 2025L, zip]),
  all(is.finite(filled$no2)), all(filled$no2 >= 0))

splice <- function(previous, current) {
  stopifnot(setequal(names(previous), names(current)))
  setorder(rbindlist(list(previous[year != 2025L], current), use.names = TRUE), zip, year, month)
}

raw_all_path <- file.path(out_dir, "all_no2_zip.csv.gz")
filled_all_path <- file.path(filled_dir, "all_no2_zip_filled.csv.gz")
raw_all <- splice(read_gz(raw_all_path), raw)
filled_all <- splice(read_gz(filled_all_path), filled)
dat <- copy(old)
target <- which(dat$year == 2025L)
aligned <- filled[match(dat$zip[target], zip)]
for (column in c("no2", "value_source", "fill_distance_km", "fill_distance_m", "fill_cell")) {
  set(dat, i = target, j = column, value = as.vector(aligned[[column]], mode = typeof(dat[[column]])))
}
stopifnot(nrow(dat) == 699300L, identical(sort(unique(dat$year)), 2005:2025),
  !anyDuplicated(dat[, .(zip, year)]), all(is.finite(dat$no2)))

stage <- tempfile(fileext = ".parquet")
arrow::write_parquet(dat, stage, compression = "zstd", compression_level = 9, use_dictionary = TRUE)
reopened <- as.data.table(arrow::read_parquet(stage))
stopifnot(identical(as.list(dat), as.list(reopened)),
  identical(as.list(historical), as.list(reopened[year != 2025L])))
if (!file.rename(stage, asset)) stop("Could not finalize ", asset)

fwrite(raw, file.path(out_dir, period))
fwrite(filled, file.path(filled_dir, period))
fwrite(raw_all, raw_all_path)
fwrite(filled_all, filled_all_path)

summary_path <- file.path(filled_dir, "fill_summary.csv")
summary <- fread(summary_path)
summary[file == period, `:=`(rows = nrow(filled), zonal_mean = sum(filled$value_source == "zonal_mean"),
  nearest_zip = sum(filled$value_source == "nearest_zip"),
  max_fill_distance_km = if (any(is.finite(filled$fill_distance_km))) max(filled$fill_distance_km, na.rm = TRUE) else NA_real_)]
fwrite(summary, summary_path)
fwrite(summary, file.path(release_dir, basename(summary_path)))

audit_path <- file.path(out_dir, "fill_audit_summary.csv")
audit <- fread(audit_path)
no2_audit <- raw_all[, .(rows = .N, zips = uniqueN(zip), missing_values = sum(is.na(no2)),
  max_fill_km = if (any(is.finite(fill_distance_m))) max(fill_distance_m, na.rm = TRUE) / 1000 else NA_real_),
  by = .(pollutant, value_source)]
audit <- rbindlist(list(audit[pollutant != "no2"], no2_audit), use.names = TRUE)
fwrite(audit, audit_path)
fwrite(audit, file.path(release_dir, basename(audit_path)))

manifest_path <- file.path(release_dir, "air_pollution_zcta_parquet_manifest.csv")
manifest <- fread(manifest_path)
manifest[pollutant == "no2", `:=`(size_bytes = file.size(asset), size_mb = round(file.size(asset) / 1024^2, 1),
  sha256 = unname(tools::sha256sum(asset)))]
fwrite(manifest, manifest_path)
print(raw[, .N, by = value_source])
print(summary[file == period])
message("Verified identical 2005-2024 Parquet rows and corrected ", nrow(filled), " ZCTAs in 2025")
