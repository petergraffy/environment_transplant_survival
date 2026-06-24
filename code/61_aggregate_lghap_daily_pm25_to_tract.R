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

terraOptions(progress = as.integer(Sys.getenv("TRACT_PM25_TERRA_PROGRESS", "1")), memfrac = 0.75)
sf_use_s2(FALSE)

raw_dir <- Sys.getenv("TRACT_PM25_RAW_DIR", file.path("data", "raw", "lghap_pm25_daily"))
region_slug <- Sys.getenv("TRACT_PM25_REGION_SLUG", "cook_county_il")
tract_statefp <- Sys.getenv("TRACT_PM25_STATEFP", "17")
tract_countyfp <- Sys.getenv("TRACT_PM25_COUNTYFP", "031")
out_dir <- Sys.getenv("TRACT_PM25_OUT_DIR", file.path("data", "processed", paste0("lghap_pm25_tract_daily_", region_slug)))
cache_dir <- file.path("data", "cache")
extract_dir <- file.path(tempdir(), "lghap_pm25_tract_extract")
skip_existing <- identical(tolower(Sys.getenv("TRACT_PM25_SKIP_EXISTING_MONTH")), "true")
keep_nc <- identical(tolower(Sys.getenv("TRACT_PM25_KEEP_NC")), "true")
search_radius <- as.numeric(Sys.getenv("TRACT_PM25_FILL_SEARCH_RADIUS_M", "2000000"))
bbox_buffer <- as.numeric(Sys.getenv("TRACT_PM25_BBOX_BUFFER_DEG", "0.25"))
raw_crop_extent <- NULL
output_extent <- NULL

parse_years <- function(value, default = 2005:2021) {
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

years <- parse_years(Sys.getenv("TRACT_PM25_YEARS"), 2005:2021)
months <- parse_int_set(Sys.getenv("TRACT_PM25_MONTHS"), 1:12)
if (any(months < 1L | months > 12L)) stop("TRACT_PM25_MONTHS must be 1-12.", call. = FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

tract_url <- "https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_us_tract_500k.zip"
tract_zip <- file.path(cache_dir, basename(tract_url))
tract_dir <- file.path(cache_dir, "cb_2020_us_tract_500k")

download_tracts <- function() {
  if (!file.exists(file.path(tract_dir, "cb_2020_us_tract_500k.shp"))) {
    dir.create(tract_dir, recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(tract_zip)) {
      log_msg("Downloading 2020 Census tract boundaries")
      download.file(tract_url, tract_zip, mode = "wb", quiet = FALSE)
    }
    unzip(tract_zip, exdir = tract_dir)
  }
  tract_dir
}

read_tracts <- function() {
  shp_dir <- download_tracts()
  tr <- st_read(file.path(shp_dir, "cb_2020_us_tract_500k.shp"), quiet = TRUE)
  keep <- intersect(c("GEOID", "STATEFP", "COUNTYFP", "TRACTCE", "geometry"), names(tr))
  tr <- tr[, keep]
  tr <- tr[tr$STATEFP == tract_statefp & tr$COUNTYFP == tract_countyfp, ]
  if (!nrow(tr)) {
    stop("No tracts matched STATEFP=", tract_statefp, " COUNTYFP=", tract_countyfp, call. = FALSE)
  }
  names(tr)[names(tr) == "GEOID"] <- "tract_geoid"
  tr$tract_id <- seq_len(nrow(tr))
  st_make_valid(tr)
}

zip_files <- function() {
  files <- list.files(raw_dir, pattern = "^LGHAP\\.Global_PM25\\.D001\\.M[0-9]{6}\\.zip$", full.names = TRUE)
  dt <- data.table(path = normalizePath(files, winslash = "/", mustWork = TRUE))
  dt[, ym := sub(".*\\.M([0-9]{6})\\.zip$", "\\1", basename(path))]
  dt[, `:=`(year = as.integer(substr(ym, 1, 4)), month = as.integer(substr(ym, 5, 6)))]
  dt[year %in% years & month %in% months][order(year, month)]
}

archive_members <- function(path, year, month) {
  listing <- unzip(path, list = TRUE)
  dt <- as.data.table(listing)
  dt <- dt[grepl("^LGHAP\\.Global_PM25\\.D001\\.A[0-9]{8}\\.nc$", Name)]
  dt[, ymd := sub(".*\\.A([0-9]{8})\\.nc$", "\\1", Name)]
  dt[, date := as.Date(ymd, format = "%Y%m%d")]
  dt[as.integer(format(date, "%Y")) == year & as.integer(format(date, "%m")) == month][order(date)]
}

extract_member <- function(zip_path, member) {
  unlink(file.path(extract_dir, "*"), recursive = TRUE, force = TRUE)
  utils::unzip(zip_path, files = member, exdir = extract_dir)
  file.path(extract_dir, member)
}

set_region_extent <- function(tracts) {
  bbox <- st_bbox(st_transform(tracts, "EPSG:4326"))
  xmin <- max(-180, as.numeric(bbox[["xmin"]]) - bbox_buffer)
  xmax <- min(180, as.numeric(bbox[["xmax"]]) + bbox_buffer)
  ymin <- max(-90, as.numeric(bbox[["ymin"]]) - bbox_buffer)
  ymax <- min(90, as.numeric(bbox[["ymax"]]) + bbox_buffer)
  list(
    raw = ext(ymin, ymax, xmin, xmax),
    output = ext(xmin, xmax, ymin, ymax)
  )
}

read_pm25_region <- function(path) {
  nc <- nc_open(path)
  on.exit(nc_close(nc), add = TRUE)

  lat <- nc$dim$lat$vals
  lon <- nc$dim$lon$vals
  lat_idx <- which(lat >= ymin(output_extent) & lat <= ymax(output_extent))
  lon_idx <- which(lon >= xmin(output_extent) & lon <= xmax(output_extent))
  vals <- ncvar_get(
    nc,
    "PM25",
    start = c(min(lat_idx), min(lon_idx)),
    count = c(length(lat_idx), length(lon_idx))
  )
  vals[vals <= -900 | vals >= 32767] <- NA_real_

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

make_zone_raster <- function(template, tracts, grid_name) {
  cache_path <- file.path(cache_dir, paste0("tract_zone_raster_", grid_name, "_v1.tif"))
  lookup_path <- file.path(cache_dir, paste0("tract_lookup_", grid_name, "_v1.csv"))
  points_path <- file.path(cache_dir, paste0("tract_points_", grid_name, "_v1.csv"))
  if (file.exists(cache_path) && file.exists(lookup_path) && file.exists(points_path)) {
    return(list(
      zones = rast(cache_path),
      lookup = fread(lookup_path, colClasses = c(tract_geoid = "character")),
      points = fread(points_path, colClasses = c(tract_geoid = "character"))
    ))
  }

  log_msg("Building tract zone raster for ", grid_name)
  if (is.na(crs(template))) crs(template) <- "EPSG:4326"
  tr_all <- st_transform(tracts, crs(template))
  tr_raster <- suppressWarnings(st_crop(tr_all, st_as_sfc(st_bbox(template))))
  tr_vect <- vect(tr_raster)
  zones <- rasterize(tr_vect, template, field = "tract_id", touches = TRUE)
  pts <- st_point_on_surface(tr_all)
  xy <- st_coordinates(pts)
  point_dt <- data.table(
    tract_id = pts$tract_id,
    tract_geoid = pts$tract_geoid,
    x = xy[, "X"],
    y = xy[, "Y"]
  )

  writeRaster(zones, cache_path, overwrite = TRUE, wopt = list(gdal = "COMPRESS=LZW"))
  fwrite(st_drop_geometry(tr_all[, c("tract_id", "tract_geoid")]), lookup_path)
  fwrite(point_dt, points_path)

  list(
    zones = rast(cache_path),
    lookup = fread(lookup_path, colClasses = c(tract_geoid = "character")),
    points = fread(points_path, colClasses = c(tract_geoid = "character"))
  )
}

fill_nearest_raster <- function(dt, r, zones) {
  missing <- which(is.na(dt$pm25_ug_m3))
  if (!length(missing)) return(dt)

  point_dt <- zones$points[match(dt$tract_geoid[missing], zones$points$tract_geoid)]
  point_vect <- vect(point_dt[, .(x, y)], geom = c("x", "y"), crs = crs(r))
  nearest <- as.data.table(extract(r, point_vect, ID = FALSE, search_radius = search_radius))
  value_col <- names(r)[1]
  fillable <- !is.na(nearest[[value_col]])

  if (any(!fillable)) {
    unresolved <- which(!fillable)
    r_ext <- ext(r)
    clamped_points <- copy(point_dt[unresolved])
    clamped_points[, x := pmin(pmax(x, xmin(r_ext)), xmax(r_ext))]
    clamped_points[, y := pmin(pmax(y, ymin(r_ext)), ymax(r_ext))]
    clamped_vect <- vect(clamped_points[, .(x, y)], geom = c("x", "y"), crs = crs(r))
    clamped_nearest <- as.data.table(extract(r, clamped_vect, ID = FALSE, search_radius = search_radius))
    clamped_fillable <- !is.na(clamped_nearest[[value_col]])
    if (any(clamped_fillable)) {
      nearest[unresolved[clamped_fillable], (value_col) := clamped_nearest[[value_col]][clamped_fillable]]
      nearest[unresolved[clamped_fillable], distance := clamped_nearest$distance[clamped_fillable]]
      nearest[unresolved[clamped_fillable], cell := clamped_nearest$cell[clamped_fillable]]
      fillable <- !is.na(nearest[[value_col]])
    }
  }

  idx <- missing[fillable]
  dt$pm25_ug_m3[idx] <- nearest[[value_col]][fillable]
  dt$value_source[missing] <- "nearest_raster_cell"
  dt$value_source[missing[!fillable]] <- "missing"
  dt$fill_distance_m[idx] <- nearest$distance[fillable]
  dt$fill_cell[idx] <- nearest$cell[fillable]
  dt
}

zonal_mean <- function(r, zones) {
  z <- crop(zones$zones, r)
  names(z) <- "tract_id"
  zs <- as.data.table(zonal(r, z, fun = "mean", na.rm = TRUE))
  setnames(zs, "tract_id", "tract_id")
  dt <- merge(zones$lookup, zs, by = "tract_id", all.x = TRUE, sort = FALSE)[, tract_id := NULL][]
  dt[, value_source := fifelse(is.na(pm25_ug_m3), "missing", "zonal_mean")]
  dt[, fill_distance_m := NA_real_]
  dt[, fill_cell := NA_integer_]
  fill_nearest_raster(dt, r, zones)
}

write_month <- function(month_dt, year, month) {
  month_dt[, `:=`(
    pollutant = "pm25",
    year = as.integer(year),
    month = as.integer(month),
    temporal_resolution = "daily"
  )]
  setcolorder(month_dt, c("tract_geoid", "date", "pollutant", "year", "month", "temporal_resolution"))
  out <- file.path(out_dir, sprintf("lghap_pm25_tract_daily_%04d_%02d.csv.gz", year, month))
  fwrite(month_dt, out)
  out
}

log_msg("Raw daily PM2.5 directory: ", raw_dir)
log_msg("Region: ", region_slug, " STATEFP=", tract_statefp, " COUNTYFP=", tract_countyfp)
log_msg("Years: ", paste(years, collapse = ", "))
log_msg("Months: ", paste(months, collapse = ", "))
files <- zip_files()
if (!nrow(files)) stop("No LGHAP daily PM2.5 monthly zip files found.", call. = FALSE)
log_msg("Monthly zip files in range: ", nrow(files))

tracts <- read_tracts()
region_extent <- set_region_extent(tracts)
raw_crop_extent <- region_extent$raw
output_extent <- region_extent$output
log_msg("Tracts in region: ", nrow(tracts))
log_msg("Raster crop extent: lon ", round(xmin(output_extent), 4), " to ", round(xmax(output_extent), 4),
        ", lat ", round(ymin(output_extent), 4), " to ", round(ymax(output_extent), 4))
manifest <- list()
zones <- NULL

for (i in seq_len(nrow(files))) {
  f <- files[i]
  out_path <- file.path(out_dir, sprintf("lghap_pm25_tract_daily_%04d_%02d.csv.gz", f$year, f$month))
  if (skip_existing && file.exists(out_path)) {
    log_msg("Skipping existing ", f$year, "-", sprintf("%02d", f$month))
    manifest[[length(manifest) + 1L]] <- data.table(
      year = f$year, month = f$month, file = basename(out_path), path = normalizePath(out_path, winslash = "/", mustWork = TRUE)
    )
    next
  }

  members <- archive_members(f$path, f$year, f$month)
  log_msg("Processing ", f$year, "-", sprintf("%02d", f$month), " days=", nrow(members), " (", i, "/", nrow(files), ")")
  day_list <- vector("list", nrow(members))
  for (j in seq_len(nrow(members))) {
    m <- members[j]
    log_msg("  day ", j, "/", nrow(members), " ", as.character(m$date))
    nc_path <- extract_member(f$path, m$Name)
    r <- read_pm25_region(nc_path)
    if (is.null(zones)) zones <- make_zone_raster(r, tracts, paste0("lghap_pm25_001deg_ncdf_", region_slug))
    dt <- zonal_mean(r, zones)
    dt[, date := m$date]
    day_list[[j]] <- dt
    if (!keep_nc) unlink(nc_path, force = TRUE)
    rm(dt, r)
    gc()
  }

  month_dt <- rbindlist(day_list, fill = TRUE)
  out <- write_month(month_dt, f$year, f$month)
  qc <- month_dt[, .(
    rows = .N,
    tracts = uniqueN(tract_geoid),
    dates = uniqueN(date),
    missing_values = sum(is.na(pm25_ug_m3)),
    zonal_mean_values = sum(value_source == "zonal_mean", na.rm = TRUE),
    nearest_fill_values = sum(value_source == "nearest_raster_cell", na.rm = TRUE),
    mean = mean(pm25_ug_m3, na.rm = TRUE),
    p01 = as.numeric(quantile(pm25_ug_m3, 0.01, na.rm = TRUE)),
    p50 = as.numeric(quantile(pm25_ug_m3, 0.50, na.rm = TRUE)),
    p99 = as.numeric(quantile(pm25_ug_m3, 0.99, na.rm = TRUE))
  )]
  qc[, `:=`(year = f$year, month = f$month, output = normalizePath(out, winslash = "/", mustWork = TRUE))]
  setcolorder(qc, c("year", "month"))
  fwrite(qc, file.path(out_dir, sprintf("lghap_pm25_tract_daily_qc_%04d_%02d.csv", f$year, f$month)))

  manifest[[length(manifest) + 1L]] <- data.table(
    year = f$year,
    month = f$month,
    file = basename(out),
    path = normalizePath(out, winslash = "/", mustWork = TRUE),
    rows = nrow(month_dt),
    tracts = uniqueN(month_dt$tract_geoid),
    dates = uniqueN(month_dt$date),
    missing_values = sum(is.na(month_dt$pm25_ug_m3))
  )
  rm(day_list, month_dt)
  gc()
}

manifest_dt <- rbindlist(manifest, fill = TRUE)
setorder(manifest_dt, year, month)
fwrite(manifest_dt, file.path(out_dir, "lghap_pm25_tract_daily_manifest.csv"))

qc_files <- list.files(out_dir, pattern = "^lghap_pm25_tract_daily_qc_[0-9]{4}_[0-9]{2}\\.csv$", full.names = TRUE)
if (length(qc_files)) {
  qc_all <- rbindlist(lapply(qc_files, fread), fill = TRUE)
  setorder(qc_all, year, month)
  fwrite(qc_all, file.path(out_dir, "lghap_pm25_tract_daily_qc_all_months.csv"))
}

log_msg("Done. Outputs written to ", normalizePath(out_dir, winslash = "/"))
