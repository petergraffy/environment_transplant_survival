user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", "4.5")
if (dir.exists(user_lib)) {
  .libPaths(unique(c(user_lib, .libPaths())))
}

library(data.table)
library(sf)

sf_use_s2(FALSE)

out_dir <- "output/zip_pollution"
filled_dir <- file.path(out_dir, "filled")
dir.create(filled_dir, recursive = TRUE, showWarnings = FALSE)

zcta_path <- "data/cache/cb_2020_us_zcta520_500k/cb_2020_us_zcta520_500k.shp"
zcta <- st_read(zcta_path, quiet = TRUE)[, c("ZCTA5CE20", "geometry")]
names(zcta)[1] <- "zip"
conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326))
conus_bbox <- st_transform(conus_bbox, st_crs(zcta))
zcta <- suppressWarnings(st_crop(zcta, conus_bbox))
zcta <- st_transform(zcta, 5070)
zcta$zip <- as.character(zcta$zip)
zcta <- zcta[order(zcta$zip), ]

centroids <- suppressWarnings(st_point_on_surface(zcta))
zips <- zcta$zip

read_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  as.data.table(read.csv(con, check.names = FALSE, colClasses = c(zip = "character")))
}

fill_one <- function(dt) {
  audit_cols <- c("value_source", "fill_distance_m", "fill_cell", "fill_distance_km")
  value_col <- setdiff(names(dt), c("zip", "pollutant", "year", "month", "temporal_resolution", audit_cols))
  if (length(value_col) != 1L) {
    stop("Expected one value column; found: ", paste(value_col, collapse = ", "))
  }

  period <- unique(dt[, .(pollutant, year, month, temporal_resolution)])
  if (nrow(period) != 1L) {
    stop("Expected one pollutant-period per file.")
  }

  full <- merge(data.table(zip = zips), dt, by = "zip", all.x = TRUE, sort = FALSE)
  full[, `:=`(
    pollutant = period$pollutant,
    year = period$year,
    month = period$month,
    temporal_resolution = period$temporal_resolution
  )]

  missing_idx <- which(is.na(full[[value_col]]))
  valid_idx <- which(!is.na(full[[value_col]]))
  full[, value_source := "zonal_mean"]
  full[, fill_distance_km := NA_real_]

  if (length(missing_idx) > 0L) {
    nearest_valid_pos <- st_nearest_feature(centroids[missing_idx, ], centroids[valid_idx, ])
    donor_idx <- valid_idx[nearest_valid_pos]
    full[[value_col]][missing_idx] <- full[[value_col]][donor_idx]
    full$value_source[missing_idx] <- "nearest_zip"
    full$fill_distance_km[missing_idx] <- as.numeric(st_distance(
      centroids[missing_idx, ],
      centroids[donor_idx, ],
      by_element = TRUE
    )) / 1000
  }

  setcolorder(full, c(
    "zip", "pollutant", "year", "month", "temporal_resolution",
    value_col, "value_source", "fill_distance_km"
  ))
  full[order(zip)]
}

period_files <- list.files(out_dir, pattern = "^(pm25|o3|no2)_[0-9]{4}(_[0-9]{2})?\\.csv\\.gz$", full.names = TRUE)
summary_rows <- vector("list", length(period_files))

for (i in seq_along(period_files)) {
  src <- period_files[[i]]
  dt <- read_gz(src)
  filled <- fill_one(dt)
  dest <- file.path(filled_dir, basename(src))
  fwrite(filled, dest)
  summary_rows[[i]] <- filled[, .(
    file = basename(src),
    rows = .N,
    zonal_mean = sum(value_source == "zonal_mean"),
    nearest_zip = sum(value_source == "nearest_zip"),
    max_fill_distance_km = suppressWarnings(max(fill_distance_km, na.rm = TRUE))
  )]
}

fill_summary <- rbindlist(summary_rows)
fill_summary[is.infinite(max_fill_distance_km), max_fill_distance_km := NA_real_]
fwrite(fill_summary, file.path(filled_dir, "fill_summary.csv"))

combine_pollutant <- function(prefix) {
  files <- list.files(filled_dir, pattern = paste0("^", prefix, "_[0-9]{4}(_[0-9]{2})?\\.csv\\.gz$"), full.names = TRUE)
  dt <- rbindlist(lapply(files, read_gz), fill = TRUE)
  setorder(dt, zip, year, month)
  fwrite(dt, file.path(filled_dir, paste0("all_", prefix, "_zip_filled.csv.gz")))
}

combine_pollutant("pm25")
combine_pollutant("o3")
combine_pollutant("no2")

message("Filled outputs written to ", normalizePath(filled_dir, winslash = "/"))
