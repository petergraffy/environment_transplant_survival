#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(data.table)
  library(ncdf4)
  library(sf)
  library(terra)
})

terraOptions(progress = as.integer(Sys.getenv("LGHAP_DAILY_TERRA_PROGRESS", "0")), memfrac = 0.75)
sf_use_s2(FALSE)

outer_zip <- Sys.getenv("LGHAP_DAILY_PM25_ZIP", "C:/Users/Peter Graffy/Downloads/8313615.zip")
target_year <- as.integer(Sys.getenv("LGHAP_DAILY_PM25_YEAR", "2021"))
target_months <- as.integer(strsplit(Sys.getenv("LGHAP_DAILY_PM25_MONTHS", "1:12"), ":", fixed = TRUE)[[1]])
if (length(target_months) == 2L) target_months <- seq(target_months[1], target_months[2])

raw_dir <- file.path("data", "raw", "lghap_pm25_daily")
out_dir <- file.path("data", "processed", "lghap_pm25_zcta_daily")
cache_dir <- file.path("data", "cache")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

skip_existing_month <- identical(tolower(Sys.getenv("LGHAP_DAILY_PM25_SKIP_EXISTING_MONTH", "false")), "true")
keep_daily_nc <- identical(tolower(Sys.getenv("LGHAP_DAILY_PM25_KEEP_NC", "false")), "true")
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

month_member <- function(year, month) {
  sprintf("LGHAP.Global_PM25.D001.M%04d%02d.zip", year, month)
}

day_member <- function(date) {
  sprintf("LGHAP.Global_PM25.D001.A%s.nc", format(date, "%Y%m%d"))
}

month_zip_path <- function(year, month) {
  file.path(raw_dir, month_member(year, month))
}

ensure_month_zip <- function(year, month) {
  dest <- month_zip_path(year, month)
  if (file.exists(dest) && file.info(dest)$size > 0) return(dest)
  if (!file.exists(outer_zip)) stop("Daily PM2.5 archive not found: ", outer_zip, call. = FALSE)
  log_msg("Extracting ", month_member(year, month), " from ", outer_zip)
  unzip(outer_zip, files = month_member(year, month), exdir = raw_dir, overwrite = TRUE)
  if (!file.exists(dest)) stop("Expected extracted monthly archive not found: ", dest, call. = FALSE)
  dest
}

extract_day_nc <- function(month_zip, date) {
  member <- day_member(date)
  dest <- file.path(raw_dir, member)
  if (file.exists(dest) && file.info(dest)$size > 0) return(dest)
  unzip(month_zip, files = member, exdir = raw_dir, overwrite = TRUE)
  if (!file.exists(dest)) stop("Expected extracted daily file not found: ", dest, call. = FALSE)
  dest
}

read_lghap_conus <- function(path) {
  nc <- nc_open(path)
  on.exit(nc_close(nc), add = TRUE)

  lat <- nc$dim$lat$vals
  lon <- nc$dim$lon$vals
  lat_idx <- which(lat >= conus_bbox["ymin"] & lat <= conus_bbox["ymax"])
  lon_idx <- which(lon >= conus_bbox["xmin"] & lon <= conus_bbox["xmax"])
  vals <- ncvar_get(
    nc,
    "PM25",
    start = c(min(lat_idx), min(lon_idx)),
    count = c(length(lat_idx), length(lon_idx))
  )
  vals[vals <= -900] <- NA_real_

  lat_sel <- lat[lat_idx]
  lon_sel <- lon[lon_idx]
  xres <- median(diff(lon_sel))
  yres <- median(diff(lat_sel))
  r <- rast(
    nrows = length(lat_sel),
    ncols = length(lon_sel),
    xmin = min(lon_sel) - xres / 2,
    xmax = max(lon_sel) + xres / 2,
    ymin = min(lat_sel) - yres / 2,
    ymax = max(lat_sel) + yres / 2,
    crs = "EPSG:4326"
  )
  values(r) <- as.vector(t(vals[rev(seq_len(nrow(vals))), , drop = FALSE]))
  names(r) <- "pm25_ug_m3"
  r
}

make_zone_raster <- function(template, zcta) {
  cache_path <- file.path(cache_dir, "zcta_zone_raster_lghap_pm25_001deg_v1.tif")
  lookup_path <- file.path(cache_dir, "zcta_lookup_lghap_pm25_001deg_v1.csv")
  points_path <- file.path(cache_dir, "zcta_points_lghap_pm25_001deg_v1.csv")
  if (file.exists(cache_path) && file.exists(lookup_path) && file.exists(points_path)) {
    return(list(
      zones = rast(cache_path),
      lookup = fread(lookup_path, colClasses = c(zip = "character")),
      points = fread(points_path, colClasses = c(zip = "character"))
    ))
  }

  log_msg("Building ZCTA zone raster for LGHAP PM2.5 0.01-degree grid")
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
  missing <- which(is.na(dt$pm25_ug_m3))
  dt[, `:=`(
    value_source = fifelse(is.na(pm25_ug_m3), "missing", "zonal_mean"),
    fill_distance_m = NA_real_,
    fill_cell = NA_integer_
  )]
  if (!length(missing)) return(dt)

  point_dt <- zones$points[match(dt$zip[missing], zones$points$zip)]
  point_vect <- vect(point_dt[, .(x, y)], geom = c("x", "y"), crs = crs(r))
  nearest <- as.data.table(extract(r, point_vect, ID = FALSE, search_radius = search_radius))
  fillable <- !is.na(nearest$pm25_ug_m3)
  if (any(fillable)) {
    idx <- missing[fillable]
    dt$pm25_ug_m3[idx] <- nearest$pm25_ug_m3[fillable]
    dt$value_source[idx] <- "nearest_raster_cell"
    dt$fill_distance_m[idx] <- nearest$distance[fillable]
    dt$fill_cell[idx] <- nearest$cell[fillable]
  }
  dt$value_source[missing[!fillable]] <- "missing"
  dt
}

aggregate_day <- function(date, month_zip, zcta) {
  nc_path <- extract_day_nc(month_zip, date)
  r <- read_lghap_conus(nc_path)
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
    pollutant = "pm25_lghap_daily",
    temporal_resolution = "daily"
  )]
  setcolorder(dt, c("zip", "date", "pollutant", "year", "month", "temporal_resolution", "pm25_ug_m3", "value_source", "fill_distance_m", "fill_cell"))
  if (!keep_daily_nc && file.exists(nc_path)) file.remove(nc_path)
  dt[]
}

process_month <- function(month, zcta) {
  dest <- file.path(out_dir, sprintf("lghap_pm25_zcta_daily_%04d_%02d.csv.gz", target_year, month))
  if (skip_existing_month && file.exists(dest)) {
    log_msg("Skipping existing month output ", target_year, "-", sprintf("%02d", month))
    return(data.table(year = target_year, month = month, output = normalizePath(dest, winslash = "/", mustWork = TRUE), skipped = TRUE))
  }

  month_zip <- ensure_month_zip(target_year, month)
  dates <- seq.Date(
    as.Date(sprintf("%04d-%02d-01", target_year, month)),
    by = "day",
    length.out = as.integer(format(seq(as.Date(sprintf("%04d-%02d-01", target_year, month)), by = "month", length.out = 2)[2] - 1, "%d"))
  )

  log_msg("Processing ", target_year, "-", sprintf("%02d", month), " days=", length(dates))
  month_dt <- rbindlist(lapply(seq_along(dates), function(i) {
    log_msg("  day ", i, "/", length(dates), " ", dates[i])
    aggregate_day(dates[i], month_zip, zcta)
  }), fill = TRUE)
  setorder(month_dt, zip, date)
  fwrite(month_dt, dest)

  qc <- data.table(
    year = target_year,
    month = month,
    rows = nrow(month_dt),
    zips = uniqueN(month_dt$zip),
    dates = uniqueN(month_dt$date),
    missing_values = sum(is.na(month_dt$pm25_ug_m3)),
    zonal_mean_values = sum(month_dt$value_source == "zonal_mean", na.rm = TRUE),
    nearest_fill_values = sum(month_dt$value_source == "nearest_raster_cell", na.rm = TRUE),
    mean = mean(month_dt$pm25_ug_m3, na.rm = TRUE),
    p01 = quantile(month_dt$pm25_ug_m3, 0.01, na.rm = TRUE),
    p50 = quantile(month_dt$pm25_ug_m3, 0.50, na.rm = TRUE),
    p99 = quantile(month_dt$pm25_ug_m3, 0.99, na.rm = TRUE),
    output = normalizePath(dest, winslash = "/", mustWork = TRUE)
  )
  fwrite(qc, file.path(out_dir, sprintf("lghap_pm25_zcta_daily_qc_%04d_%02d.csv", target_year, month)))
  log_msg("Wrote ", dest, " rows=", nrow(month_dt), " missing=", qc$missing_values)
  qc[, .(year, month, output, skipped = FALSE)]
}

log_msg("LGHAP daily PM2.5 archive: ", outer_zip)
log_msg("LGHAP daily PM2.5 year: ", target_year)
log_msg("Months: ", paste(target_months, collapse = ", "))
zcta <- read_zctas()
manifest <- rbindlist(lapply(target_months, process_month, zcta = zcta), fill = TRUE)
fwrite(manifest, file.path(out_dir, sprintf("lghap_pm25_zcta_daily_manifest_%04d.csv", target_year)))
qc_files <- list.files(out_dir, pattern = sprintf("^lghap_pm25_zcta_daily_qc_%04d_[0-9]{2}\\.csv$", target_year), full.names = TRUE)
if (length(qc_files)) {
  qc_all <- rbindlist(lapply(qc_files, fread), fill = TRUE)
  setorder(qc_all, year, month)
  fwrite(qc_all, file.path(out_dir, sprintf("lghap_pm25_zcta_daily_qc_%04d_all_months.csv", target_year)))
}
log_msg("Done. Outputs written to ", normalizePath(out_dir, winslash = "/"))
