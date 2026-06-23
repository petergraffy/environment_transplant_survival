#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(terra)
})

terraOptions(progress = as.integer(Sys.getenv("DAILY_O3_TERRA_PROGRESS", "0")), memfrac = 0.75)
sf_use_s2(FALSE)

target_year <- as.integer(Sys.getenv("DAILY_O3_YEAR", "2024"))
target_months <- as.integer(strsplit(Sys.getenv("DAILY_O3_MONTHS", "1:12"), ":", fixed = TRUE)[[1]])
if (length(target_months) == 2L) target_months <- seq(target_months[1], target_months[2])

outer_zip <- Sys.getenv("DAILY_O3_ZIP", sprintf("C:/Users/Peter Graffy/Downloads/%04d.zip", target_year))
raw_dir <- file.path("data", "raw", "daily_o3")
out_dir <- file.path("data", "processed", "o3_zcta_daily")
cache_dir <- file.path("data", "cache")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

skip_existing_month <- identical(tolower(Sys.getenv("DAILY_O3_SKIP_EXISTING_MONTH", "false")), "true")
keep_daily_nc <- identical(tolower(Sys.getenv("DAILY_O3_KEEP_NC", "false")), "true")
conus_bbox <- c(xmin = -125, xmax = -66, ymin = 24, ymax = 50)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

zcta_url <- "https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_us_zcta520_500k.zip"
zcta_zip <- file.path(cache_dir, basename(zcta_url))
zcta_dir <- file.path(cache_dir, "cb_2020_us_zcta520_500k")

download_zctas <- function() {
  if (!file.exists(file.path(zcta_dir, "cb_2020_us_zcta520_500k.shp"))) {
    dir.create(zcta_dir, recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(zcta_zip)) {
      log_msg("Downloading 2020 Census ZCTA boundaries")
      download.file(zcta_url, zcta_zip, mode = "wb", quiet = FALSE)
    }
    unzip(zcta_zip, exdir = zcta_dir)
  }
  zcta_dir
}

read_zctas <- function() {
  shp_dir <- download_zctas()
  z <- st_read(file.path(shp_dir, "cb_2020_us_zcta520_500k.shp"), quiet = TRUE)
  z <- z[, c("ZCTA5CE20", "geometry")]
  names(z)[1] <- "zip"
  z$zip_id <- seq_len(nrow(z))
  st_make_valid(z)
}

day_member <- function(date) {
  sprintf("%04d/%s.nc", as.integer(format(date, "%Y")), format(date, "%Y%m%d"))
}

day_nc_path <- function(date) {
  file.path(raw_dir, basename(day_member(date)))
}

extract_day_nc <- function(date) {
  dest <- day_nc_path(date)
  if (file.exists(dest) && file.info(dest)$size > 0) return(dest)
  if (!file.exists(outer_zip)) stop("Daily ozone archive not found: ", outer_zip, call. = FALSE)

  member <- day_member(date)
  members <- unzip(outer_zip, list = TRUE)$Name
  hit <- members[members == member]
  if (!length(hit)) hit <- members[basename(members) == basename(member)]
  if (!length(hit)) stop("Could not find daily ozone NetCDF in archive for ", date, call. = FALSE)

  unzip(outer_zip, files = hit[1], exdir = raw_dir, overwrite = TRUE)
  extracted <- file.path(raw_dir, hit[1])
  if (!file.exists(extracted)) stop("Expected extracted daily ozone file not found: ", extracted, call. = FALSE)
  if (normalizePath(extracted, winslash = "/", mustWork = TRUE) != normalizePath(dest, winslash = "/", mustWork = FALSE)) {
    file.copy(extracted, dest, overwrite = TRUE)
  }
  dest
}

read_daily_o3 <- function(path) {
  r <- rast(path)
  if (!"mda8" %in% names(r)) {
    idx <- grep("mda8|o3|ozone", names(r), ignore.case = TRUE)
    if (!length(idx)) stop("No ozone/MDA8 layer found in ", path, call. = FALSE)
    r <- r[[idx[1]]]
  } else {
    r <- r[["mda8"]]
  }
  if (is.na(crs(r))) crs(r) <- "EPSG:4326"
  r <- crop(r, ext(conus_bbox["xmin"], conus_bbox["xmax"], conus_bbox["ymin"], conus_bbox["ymax"]))
  r[r <= -900] <- NA
  names(r) <- "o3_ppb"
  r
}

make_zone_raster <- function(template, zcta) {
  cache_path <- file.path(cache_dir, "zcta_zone_raster_daily_o3_001deg_v1.tif")
  lookup_path <- file.path(cache_dir, "zcta_lookup_daily_o3_001deg_v1.csv")
  points_path <- file.path(cache_dir, "zcta_points_daily_o3_001deg_v1.csv")
  if (file.exists(cache_path) && file.exists(lookup_path) && file.exists(points_path)) {
    return(list(
      zones = rast(cache_path),
      lookup = fread(lookup_path, colClasses = c(zip = "character")),
      points = fread(points_path, colClasses = c(zip = "character"))
    ))
  }

  log_msg("Building ZCTA zone raster for daily ozone 0.01-degree grid")
  template <- template[[1]]
  z_all <- st_transform(zcta, crs(template))
  bbox <- st_as_sfc(st_bbox(c(
    xmin = unname(conus_bbox["xmin"]),
    ymin = unname(conus_bbox["ymin"]),
    xmax = unname(conus_bbox["xmax"]),
    ymax = unname(conus_bbox["ymax"])
  ), crs = 4326))
  bbox <- st_transform(bbox, st_crs(z_all))
  z_all <- suppressWarnings(st_crop(z_all, bbox))
  z_raster <- suppressWarnings(st_crop(z_all, st_as_sfc(st_bbox(template))))
  zones <- rasterize(vect(z_raster), template, field = "zip_id", touches = TRUE)
  pts <- suppressWarnings(st_point_on_surface(z_all))
  xy <- st_coordinates(pts)
  point_dt <- data.table(zip_id = pts$zip_id, zip = pts$zip, x = xy[, "X"], y = xy[, "Y"])
  writeRaster(zones, cache_path, overwrite = TRUE, wopt = list(gdal = "COMPRESS=LZW"))
  fwrite(st_drop_geometry(z_all[, c("zip_id", "zip")]), lookup_path)
  fwrite(point_dt, points_path)
  list(
    zones = rast(cache_path),
    lookup = fread(lookup_path, colClasses = c(zip = "character")),
    points = fread(points_path, colClasses = c(zip = "character"))
  )
}

fill_nearest_raster <- function(dt, r, zones, search_radius = 500000) {
  missing <- which(is.na(dt$o3_ppb))
  dt[, `:=`(
    value_source = fifelse(is.na(o3_ppb), "missing", "zonal_mean"),
    fill_distance_m = NA_real_,
    fill_cell = NA_integer_
  )]
  if (!length(missing)) return(dt)

  point_dt <- zones$points[match(dt$zip[missing], zones$points$zip)]
  point_vect <- vect(point_dt[, .(x, y)], geom = c("x", "y"), crs = crs(r))
  nearest <- as.data.table(extract(r, point_vect, ID = FALSE, search_radius = search_radius))
  fillable <- !is.na(nearest$o3_ppb)
  if (any(fillable)) {
    idx <- missing[fillable]
    dt$o3_ppb[idx] <- nearest$o3_ppb[fillable]
    dt$value_source[idx] <- "nearest_raster_cell"
    dt$fill_distance_m[idx] <- nearest$distance[fillable]
    dt$fill_cell[idx] <- nearest$cell[fillable]
  }
  dt$value_source[missing[!fillable]] <- "missing"
  dt
}

aggregate_day <- function(date, zcta) {
  nc_path <- extract_day_nc(date)
  r <- read_daily_o3(nc_path)
  zones <- make_zone_raster(r, zcta)
  z <- crop(zones$zones, r)
  names(z) <- "zip_id"
  zs <- as.data.table(zonal(r, z, fun = "mean", na.rm = TRUE))
  dt <- merge(zones$lookup, zs, by = "zip_id", all.x = TRUE, sort = FALSE)[, zip_id := NULL][]
  dt <- fill_nearest_raster(dt, r, zones)
  dt[, `:=`(
    date = date,
    year = as.integer(format(date, "%Y")),
    month = as.integer(format(date, "%m")),
    pollutant = "o3_daily_mda8",
    temporal_resolution = "daily"
  )]
  setcolorder(dt, c("zip", "date", "pollutant", "year", "month", "temporal_resolution", "o3_ppb", "value_source", "fill_distance_m", "fill_cell"))
  if (!keep_daily_nc && file.exists(nc_path)) file.remove(nc_path)
  nested_dir <- file.path(raw_dir, as.character(target_year))
  if (!keep_daily_nc && dir.exists(nested_dir)) unlink(nested_dir, recursive = TRUE, force = TRUE)
  dt[]
}

month_dates <- function(year, month) {
  start <- as.Date(sprintf("%04d-%02d-01", year, month))
  end <- seq(start, by = "month", length.out = 2)[2] - 1
  seq.Date(start, end, by = "day")
}

process_month <- function(month, zcta) {
  dest <- file.path(out_dir, sprintf("o3_zcta_daily_%04d_%02d.csv.gz", target_year, month))
  if (skip_existing_month && file.exists(dest)) {
    log_msg("Skipping existing month output ", target_year, "-", sprintf("%02d", month))
    return(data.table(year = target_year, month = month, output = normalizePath(dest, winslash = "/", mustWork = TRUE), skipped = TRUE))
  }

  dates <- month_dates(target_year, month)
  log_msg("Processing daily ozone ", target_year, "-", sprintf("%02d", month), " days=", length(dates))
  month_dt <- rbindlist(lapply(seq_along(dates), function(i) {
    log_msg("  day ", i, "/", length(dates), " ", dates[i])
    aggregate_day(dates[i], zcta)
  }), fill = TRUE)
  setorder(month_dt, zip, date)
  fwrite(month_dt, dest)

  qc <- data.table(
    year = target_year,
    month = month,
    rows = nrow(month_dt),
    zips = uniqueN(month_dt$zip),
    dates = uniqueN(month_dt$date),
    missing_values = sum(is.na(month_dt$o3_ppb)),
    zonal_mean_values = sum(month_dt$value_source == "zonal_mean", na.rm = TRUE),
    nearest_fill_values = sum(month_dt$value_source == "nearest_raster_cell", na.rm = TRUE),
    mean = mean(month_dt$o3_ppb, na.rm = TRUE),
    p01 = quantile(month_dt$o3_ppb, 0.01, na.rm = TRUE),
    p50 = quantile(month_dt$o3_ppb, 0.50, na.rm = TRUE),
    p99 = quantile(month_dt$o3_ppb, 0.99, na.rm = TRUE),
    output = normalizePath(dest, winslash = "/", mustWork = TRUE)
  )
  fwrite(qc, file.path(out_dir, sprintf("o3_zcta_daily_qc_%04d_%02d.csv", target_year, month)))
  log_msg("Wrote ", dest, " rows=", nrow(month_dt), " missing=", qc$missing_values)
  qc[, .(year, month, output, skipped = FALSE)]
}

log_msg("Daily ozone archive: ", outer_zip)
log_msg("Daily ozone year: ", target_year)
log_msg("Months: ", paste(target_months, collapse = ", "))
zcta <- read_zctas()
manifest <- rbindlist(lapply(target_months, process_month, zcta = zcta), fill = TRUE)
fwrite(manifest, file.path(out_dir, sprintf("o3_zcta_daily_manifest_%04d.csv", target_year)))
qc_files <- list.files(out_dir, pattern = sprintf("^o3_zcta_daily_qc_%04d_[0-9]{2}\\.csv$", target_year), full.names = TRUE)
if (length(qc_files)) {
  qc_all <- rbindlist(lapply(qc_files, fread), fill = TRUE)
  setorder(qc_all, year, month)
  fwrite(qc_all, file.path(out_dir, sprintf("o3_zcta_daily_qc_%04d_all_months.csv", target_year)))
}
log_msg("Done. Outputs written to ", normalizePath(out_dir, winslash = "/"))
