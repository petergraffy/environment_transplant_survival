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

terraOptions(progress = as.integer(Sys.getenv("TRACT_O3_TERRA_PROGRESS", "0")), memfrac = 0.75)
sf_use_s2(FALSE)

region_slug <- Sys.getenv("TRACT_O3_REGION_SLUG", "cook_county_il")
tract_statefp <- Sys.getenv("TRACT_O3_STATEFP", "17")
tract_countyfp <- Sys.getenv("TRACT_O3_COUNTYFP", "031")
download_dir <- Sys.getenv("TRACT_O3_DOWNLOAD_DIR", "C:/Users/Peter Graffy/Downloads")
raw_dir <- Sys.getenv("TRACT_O3_RAW_DIR", file.path("data", "raw", "daily_o3"))
out_dir <- Sys.getenv("TRACT_O3_OUT_DIR", file.path("data", "processed", paste0("o3_tract_daily_", region_slug)))
cache_dir <- file.path("data", "cache")
extract_dir <- file.path(tempdir(), "o3_tract_extract")
skip_existing <- identical(tolower(Sys.getenv("TRACT_O3_SKIP_EXISTING_MONTH")), "true")
keep_nc <- identical(tolower(Sys.getenv("TRACT_O3_KEEP_NC")), "true")
search_radius <- as.numeric(Sys.getenv("TRACT_O3_FILL_SEARCH_RADIUS_M", "500000"))
bbox_buffer <- as.numeric(Sys.getenv("TRACT_O3_BBOX_BUFFER_DEG", "0.25"))
output_extent <- NULL

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

years <- parse_int_set(Sys.getenv("TRACT_O3_YEARS"), 2005:2024)
months <- parse_int_set(Sys.getenv("TRACT_O3_MONTHS"), 1:12)
if (any(months < 1L | months > 12L)) stop("TRACT_O3_MONTHS must be 1-12.", call. = FALSE)

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
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

set_region_extent <- function(tracts) {
  bbox <- st_bbox(st_transform(tracts, "EPSG:4326"))
  ext(
    max(-180, as.numeric(bbox[["xmin"]]) - bbox_buffer),
    min(180, as.numeric(bbox[["xmax"]]) + bbox_buffer),
    max(-90, as.numeric(bbox[["ymin"]]) - bbox_buffer),
    min(90, as.numeric(bbox[["ymax"]]) + bbox_buffer)
  )
}

zip_path <- function(year) {
  file.path(download_dir, sprintf("%04d.zip", year))
}

day_member <- function(date) {
  sprintf("%04d/%s.nc", as.integer(format(date, "%Y")), format(date, "%Y%m%d"))
}

extract_day_nc <- function(date) {
  zip_file <- zip_path(as.integer(format(date, "%Y")))
  if (!file.exists(zip_file)) stop("Daily ozone archive not found: ", zip_file, call. = FALSE)
  member <- day_member(date)
  dest <- file.path(extract_dir, basename(member))
  if (file.exists(dest) && file.info(dest)$size > 0) return(dest)

  members <- unzip(zip_file, list = TRUE)$Name
  hit <- members[members == member]
  if (!length(hit)) hit <- members[basename(members) == basename(member)]
  if (!length(hit)) stop("Could not find daily ozone NetCDF in archive for ", date, call. = FALSE)

  unlink(file.path(extract_dir, "*"), recursive = TRUE, force = TRUE)
  unzip(zip_file, files = hit[1], exdir = extract_dir, overwrite = TRUE)
  extracted <- file.path(extract_dir, hit[1])
  if (!file.exists(extracted)) stop("Expected extracted daily ozone file not found: ", extracted, call. = FALSE)
  if (normalizePath(extracted, winslash = "/", mustWork = TRUE) != normalizePath(dest, winslash = "/", mustWork = FALSE)) {
    file.copy(extracted, dest, overwrite = TRUE)
  }
  dest
}

read_daily_o3_region <- function(path) {
  r <- rast(path)
  if (!"mda8" %in% names(r)) {
    idx <- grep("mda8|o3|ozone", names(r), ignore.case = TRUE)
    if (!length(idx)) stop("No ozone/MDA8 layer found in ", path, call. = FALSE)
    r <- r[[idx[1]]]
  } else {
    r <- r[["mda8"]]
  }
  if (is.na(crs(r))) crs(r) <- "EPSG:4326"
  r <- crop(r, output_extent)
  r[r <= -900] <- NA
  names(r) <- "o3_ppb"
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
  tr_all <- st_transform(tracts, crs(template))
  tr_raster <- suppressWarnings(st_crop(tr_all, st_as_sfc(st_bbox(template))))
  zones <- rasterize(vect(tr_raster), template, field = "tract_id", touches = TRUE)
  pts <- suppressWarnings(st_point_on_surface(tr_all))
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
  missing <- which(is.na(dt$o3_ppb))
  if (!length(missing)) return(dt)
  point_dt <- zones$points[match(dt$tract_geoid[missing], zones$points$tract_geoid)]
  point_vect <- vect(point_dt[, .(x, y)], geom = c("x", "y"), crs = crs(r))
  nearest <- as.data.table(extract(r, point_vect, ID = FALSE, cells = TRUE))
  fillable <- !is.na(nearest$o3_ppb)
  if (any(fillable)) {
    idx <- missing[fillable]
    dt$o3_ppb[idx] <- nearest$o3_ppb[fillable]
    dt$value_source[idx] <- "nearest_raster_cell"
    dt$fill_distance_m[idx] <- NA_real_
    dt$fill_cell[idx] <- nearest$cell[fillable]
  }
  unresolved <- missing[!fillable]
  if (length(unresolved)) {
    grid <- as.data.table(as.data.frame(r, xy = TRUE, cells = TRUE, na.rm = TRUE))
    setnames(grid, "o3_ppb", "value")
    unresolved_points <- point_dt[!fillable]
    for (k in seq_len(nrow(unresolved_points))) {
      dx <- grid$x - unresolved_points$x[k]
      dy <- grid$y - unresolved_points$y[k]
      nearest_idx <- which.min(dx * dx + dy * dy)
      row_idx <- unresolved[k]
      dt$o3_ppb[row_idx] <- grid$value[nearest_idx]
      dt$value_source[row_idx] <- "nearest_nonmissing_raster_cell"
      dt$fill_distance_m[row_idx] <- sqrt(dx[nearest_idx] * dx[nearest_idx] + dy[nearest_idx] * dy[nearest_idx]) * 111320
      dt$fill_cell[row_idx] <- grid$cell[nearest_idx]
    }
  }
  dt$value_source[missing[is.na(dt$o3_ppb[missing])]] <- "missing"
  dt
}

zonal_mean <- function(r, zones) {
  z <- crop(zones$zones, r)
  names(z) <- "tract_id"
  zs <- as.data.table(zonal(r, z, fun = "mean", na.rm = TRUE))
  dt <- merge(zones$lookup, zs, by = "tract_id", all.x = TRUE, sort = FALSE)[, tract_id := NULL][]
  dt[, `:=`(
    value_source = fifelse(is.na(o3_ppb), "missing", "zonal_mean"),
    fill_distance_m = NA_real_,
    fill_cell = NA_integer_
  )]
  fill_nearest_raster(dt, r, zones)
}

month_dates <- function(year, month) {
  start <- as.Date(sprintf("%04d-%02d-01", year, month))
  end <- seq(start, by = "month", length.out = 2)[2] - 1
  seq.Date(start, end, by = "day")
}

write_month <- function(month_dt, year, month) {
  month_dt[, `:=`(
    pollutant = "o3_daily_mda8",
    year = as.integer(year),
    month = as.integer(month),
    temporal_resolution = "daily"
  )]
  setcolorder(month_dt, c("tract_geoid", "date", "pollutant", "year", "month", "temporal_resolution"))
  out <- file.path(out_dir, sprintf("o3_tract_daily_%04d_%02d.csv.gz", year, month))
  fwrite(month_dt, out)
  out
}

log_msg("Daily ozone archive directory: ", download_dir)
log_msg("Region: ", region_slug, " STATEFP=", tract_statefp, " COUNTYFP=", tract_countyfp)
log_msg("Years: ", paste(years, collapse = ", "))
log_msg("Months: ", paste(months, collapse = ", "))

files <- CJ(year = years, month = months)[order(year, month)]
missing_zips <- files[!file.exists(zip_path(year)), unique(zip_path(year))]
if (length(missing_zips)) stop("Missing daily ozone archives: ", paste(missing_zips, collapse = ", "), call. = FALSE)

tracts <- read_tracts()
output_extent <- set_region_extent(tracts)
log_msg("Tracts in region: ", nrow(tracts))
log_msg("Raster crop extent: lon ", round(xmin(output_extent), 4), " to ", round(xmax(output_extent), 4),
        ", lat ", round(ymin(output_extent), 4), " to ", round(ymax(output_extent), 4))

manifest <- list()
zones <- NULL

for (i in seq_len(nrow(files))) {
  f <- files[i]
  out_path <- file.path(out_dir, sprintf("o3_tract_daily_%04d_%02d.csv.gz", f$year, f$month))
  if (skip_existing && file.exists(out_path)) {
    log_msg("Skipping existing ", f$year, "-", sprintf("%02d", f$month))
    manifest[[length(manifest) + 1L]] <- data.table(
      year = f$year, month = f$month, file = basename(out_path), path = normalizePath(out_path, winslash = "/", mustWork = TRUE)
    )
    next
  }

  dates <- month_dates(f$year, f$month)
  log_msg("Processing daily ozone ", f$year, "-", sprintf("%02d", f$month), " days=", length(dates), " (", i, "/", nrow(files), ")")
  day_list <- vector("list", length(dates))
  for (j in seq_along(dates)) {
    log_msg("  day ", j, "/", length(dates), " ", dates[j])
    nc_path <- extract_day_nc(dates[j])
    r <- read_daily_o3_region(nc_path)
    if (is.null(zones)) zones <- make_zone_raster(r, tracts, paste0("daily_o3_001deg_", region_slug))
    dt <- zonal_mean(r, zones)
    dt[, date := dates[j]]
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
    missing_values = sum(is.na(o3_ppb)),
    zonal_mean_values = sum(value_source == "zonal_mean", na.rm = TRUE),
    nearest_fill_values = sum(value_source == "nearest_raster_cell", na.rm = TRUE),
    mean = mean(o3_ppb, na.rm = TRUE),
    p01 = as.numeric(quantile(o3_ppb, 0.01, na.rm = TRUE)),
    p50 = as.numeric(quantile(o3_ppb, 0.50, na.rm = TRUE)),
    p99 = as.numeric(quantile(o3_ppb, 0.99, na.rm = TRUE))
  )]
  qc[, `:=`(year = f$year, month = f$month, output = normalizePath(out, winslash = "/", mustWork = TRUE))]
  setcolorder(qc, c("year", "month"))
  fwrite(qc, file.path(out_dir, sprintf("o3_tract_daily_qc_%04d_%02d.csv", f$year, f$month)))

  manifest[[length(manifest) + 1L]] <- data.table(
    year = f$year,
    month = f$month,
    file = basename(out),
    path = normalizePath(out, winslash = "/", mustWork = TRUE),
    rows = nrow(month_dt),
    tracts = uniqueN(month_dt$tract_geoid),
    dates = uniqueN(month_dt$date),
    missing_values = sum(is.na(month_dt$o3_ppb))
  )
  rm(day_list, month_dt)
  gc()
}

manifest_dt <- rbindlist(manifest, fill = TRUE)
setorder(manifest_dt, year, month)
fwrite(manifest_dt, file.path(out_dir, "o3_tract_daily_manifest.csv"))

qc_files <- list.files(out_dir, pattern = "^o3_tract_daily_qc_[0-9]{4}_[0-9]{2}\\.csv$", full.names = TRUE)
if (length(qc_files)) {
  qc_all <- rbindlist(lapply(qc_files, fread), fill = TRUE)
  setorder(qc_all, year, month)
  fwrite(qc_all, file.path(out_dir, "o3_tract_daily_qc_all_months.csv"))
}

log_msg("Done. Outputs written to ", normalizePath(out_dir, winslash = "/"))
