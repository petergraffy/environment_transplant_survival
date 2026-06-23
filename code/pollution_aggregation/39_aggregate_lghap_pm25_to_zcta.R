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

terraOptions(progress = as.integer(Sys.getenv("LGHAP_TERRA_PROGRESS", "0")), memfrac = 0.75)
sf_use_s2(FALSE)

zip_path <- Sys.getenv("LGHAP_PM25_ZIP", "C:/Users/Peter Graffy/Downloads/11435542.zip")
years <- as.integer(strsplit(Sys.getenv("LGHAP_PM25_YEARS", "2000:2021"), ":", fixed = TRUE)[[1]])
if (length(years) == 2L) years <- seq(years[1], years[2])

raw_dir <- file.path("data", "raw", "lghap_pm25")
out_dir <- file.path("data", "processed", "lghap_pm25_zcta_annual")
cache_dir <- file.path("data", "cache")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

skip_existing_year <- identical(tolower(Sys.getenv("LGHAP_PM25_SKIP_EXISTING_YEAR", "false")), "true")
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

archive_member <- function(year) {
  sprintf("LGHAP.Global_PM25.Y001.A%04d.nc", year)
}

raw_path <- function(year) {
  file.path(raw_dir, archive_member(year))
}

ensure_raw_file <- function(year) {
  dest <- raw_path(year)
  if (file.exists(dest) && file.info(dest)$size > 0) return(dest)
  if (!file.exists(zip_path)) stop("PM2.5 archive not found: ", zip_path, call. = FALSE)
  log_msg("Extracting ", archive_member(year), " from ", zip_path)
  unzip(zip_path, files = archive_member(year), exdir = raw_dir, overwrite = TRUE)
  if (!file.exists(dest)) stop("Expected extracted file not found: ", dest, call. = FALSE)
  dest
}

read_lghap_conus <- function(path) {
  nc <- nc_open(path)
  on.exit(nc_close(nc), add = TRUE)

  lat <- nc$dim$lat$vals
  lon <- nc$dim$lon$vals
  lat_idx <- which(lat >= conus_bbox["ymin"] & lat <= conus_bbox["ymax"])
  lon_idx <- which(lon >= conus_bbox["xmin"] & lon <= conus_bbox["xmax"])
  if (!length(lat_idx) || !length(lon_idx)) stop("No CONUS cells found in ", path, call. = FALSE)

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

aggregate_year <- function(year, zcta) {
  dest <- file.path(out_dir, sprintf("lghap_pm25_zcta_annual_%04d.csv.gz", year))
  if (skip_existing_year && file.exists(dest)) {
    log_msg("Skipping existing year output ", year)
    return(data.table(year = year, output = normalizePath(dest, winslash = "/", mustWork = TRUE), skipped = TRUE))
  }

  path <- ensure_raw_file(year)
  log_msg("Reading corrected CONUS raster for ", year)
  r <- read_lghap_conus(path)
  zones <- make_zone_raster(r, zcta)
  z <- crop(zones$zones, r)
  names(z) <- "zip_id"
  log_msg("Calculating ZCTA means for ", year)
  zs <- as.data.table(zonal(r, z, fun = "mean", na.rm = TRUE))
  setnames(zs, "zip_id", "zip_id")
  dt <- merge(zones$lookup, zs, by = "zip_id", all.x = TRUE, sort = FALSE)[, zip_id := NULL][]
  dt <- fill_nearest_raster(dt, r, zones)
  dt[, `:=`(pollutant = "pm25_lghap", year = year, temporal_resolution = "annual")]
  setcolorder(dt, c("zip", "pollutant", "year", "temporal_resolution", "pm25_ug_m3", "value_source", "fill_distance_m", "fill_cell"))
  setorder(dt, zip)
  fwrite(dt, dest)
  log_msg("Wrote ", dest, " rows=", nrow(dt), " missing=", sum(is.na(dt$pm25_ug_m3)))

  qc <- data.table(
    year = year,
    rows = nrow(dt),
    zips = uniqueN(dt$zip),
    missing_values = sum(is.na(dt$pm25_ug_m3)),
    zonal_mean_values = sum(dt$value_source == "zonal_mean", na.rm = TRUE),
    nearest_fill_values = sum(dt$value_source == "nearest_raster_cell", na.rm = TRUE),
    mean = mean(dt$pm25_ug_m3, na.rm = TRUE),
    p01 = quantile(dt$pm25_ug_m3, 0.01, na.rm = TRUE),
    p50 = quantile(dt$pm25_ug_m3, 0.50, na.rm = TRUE),
    p99 = quantile(dt$pm25_ug_m3, 0.99, na.rm = TRUE),
    output = normalizePath(dest, winslash = "/", mustWork = TRUE)
  )
  fwrite(qc, file.path(out_dir, sprintf("lghap_pm25_zcta_annual_qc_%04d.csv", year)))
  qc[, .(year, output, skipped = FALSE)]
}

log_msg("LGHAP PM2.5 archive: ", zip_path)
log_msg("LGHAP PM2.5 years: ", paste(range(years), collapse = "-"))
log_msg("Note: downloaded NetCDFs contain one PM2.5 layer per year and no daily time dimension.")
zcta <- read_zctas()
manifest <- rbindlist(lapply(years, aggregate_year, zcta = zcta), fill = TRUE)
fwrite(manifest, file.path(out_dir, "lghap_pm25_zcta_annual_manifest.csv"))
qc_files <- list.files(out_dir, pattern = "^lghap_pm25_zcta_annual_qc_[0-9]{4}\\.csv$", full.names = TRUE)
if (length(qc_files)) {
  qc_all <- rbindlist(lapply(qc_files, fread), fill = TRUE)
  setorder(qc_all, year)
  fwrite(qc_all, file.path(out_dir, "lghap_pm25_zcta_annual_qc_all_years.csv"))
}
log_msg("Done. Outputs written to ", normalizePath(out_dir, winslash = "/"))
