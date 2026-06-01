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

terraOptions(progress = as.integer(Sys.getenv("GRIDMET_TERRA_PROGRESS", "0")), memfrac = 0.75)
sf_use_s2(FALSE)

gridmet_base_url <- "https://www.northwestknowledge.net/metdata/data"
years <- as.integer(strsplit(Sys.getenv("GRIDMET_YEARS", "2005:2025"), ":", fixed = TRUE)[[1]])
if (length(years) == 2L) years <- seq(years[1], years[2])

default_vars <- c("tmmx", "tmmn", "rmax", "rmin", "pr", "sph", "srad", "vs", "pet", "etr", "erc")
target_vars <- strsplit(Sys.getenv("GRIDMET_VARS", paste(default_vars, collapse = ",")), ",", fixed = TRUE)[[1]]
target_vars <- trimws(tolower(target_vars))

raw_dir <- file.path("data", "raw", "gridmet")
out_dir <- file.path("data", "processed", "gridmet_zcta_daily")
cache_dir <- file.path("data", "cache")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

skip_download <- identical(tolower(Sys.getenv("GRIDMET_SKIP_DOWNLOAD")), "true")
skip_existing_year <- identical(tolower(Sys.getenv("GRIDMET_SKIP_EXISTING_YEAR")), "true")
show_fill_progress <- identical(tolower(Sys.getenv("GRIDMET_FILL_PROGRESS", "false")), "true")

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

gridmet_url <- function(var, year) {
  sprintf("%s/%s_%04d.nc", gridmet_base_url, var, year)
}

gridmet_path <- function(var, year) {
  file.path(raw_dir, sprintf("%s_%04d.nc", var, year))
}

download_gridmet_file <- function(var, year) {
  dest <- gridmet_path(var, year)
  if (file.exists(dest) && file.info(dest)$size > 0) return(dest)
  if (skip_download) stop("Missing ", dest, " and GRIDMET_SKIP_DOWNLOAD=true", call. = FALSE)

  url <- gridmet_url(var, year)
  tmp <- paste0(dest, ".download")
  if (file.exists(tmp)) file.remove(tmp)
  log_msg("Downloading ", basename(dest))
  download.file(url, tmp, mode = "wb", quiet = FALSE)
  file.rename(tmp, dest)
  dest
}

make_zone_raster <- function(template, zcta) {
  cache_path <- file.path(cache_dir, "zcta_zone_raster_gridmet_4km_v1.tif")
  lookup_path <- file.path(cache_dir, "zcta_lookup_gridmet_4km_v1.csv")
  points_path <- file.path(cache_dir, "zcta_points_gridmet_4km_v1.csv")
  if (file.exists(cache_path) && file.exists(lookup_path) && file.exists(points_path)) {
    return(list(
      zones = rast(cache_path),
      lookup = fread(lookup_path, colClasses = c(zip = "character")),
      points = fread(points_path, colClasses = c(zip = "character"))
    ))
  }

  log_msg("Building ZCTA zone raster for gridMET 4-km grid")
  template <- template[[1]]
  if (is.na(crs(template))) crs(template) <- "EPSG:4326"
  z_all <- st_transform(zcta, crs(template))
  conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326))
  conus_bbox <- st_transform(conus_bbox, st_crs(z_all))
  z_all <- suppressWarnings(st_crop(z_all, conus_bbox))
  z_raster <- suppressWarnings(st_crop(z_all, st_as_sfc(st_bbox(template))))
  z_vect <- vect(z_raster)
  zones <- rasterize(z_vect, template, field = "zip_id", touches = TRUE)
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

layer_dates <- function(r, year) {
  n <- nlyr(r)
  dates <- as.Date(time(r))
  if (length(dates) != n || any(is.na(dates))) {
    dates <- seq.Date(as.Date(sprintf("%04d-01-01", year)), by = "day", length.out = n)
  }
  dates
}

standardize_units <- function(dt, var) {
  if (var %in% c("tmmx", "tmmn")) {
    dt[, value := value - 273.15]
  }
  dt
}

fill_nearest_by_layer <- function(dt, r, zones, var, search_radius = 500000) {
  missing <- which(is.na(dt$value))
  if (!length(missing)) {
    dt[, `:=`(value_source = "zonal_mean", fill_distance_m = NA_real_, fill_cell = NA_integer_)]
    return(dt)
  }

  dt[, `:=`(value_source = fifelse(is.na(value), "missing", "zonal_mean"), fill_distance_m = NA_real_, fill_cell = NA_integer_)]
  date_levels <- unique(dt$date[missing])
  point_cells <- copy(zones$points)
  point_cells[, cell := cellFromXY(r[[1]], as.matrix(.SD)), .SDcols = c("x", "y")]
  point_cells[, cell_distance_m := 0]
  first_vals <- values(r[[1]], mat = FALSE)
  needs_nearest <- is.na(point_cells$cell) | is.na(first_vals[point_cells$cell])
  if (any(needs_nearest)) {
    valid_cells <- which(!is.na(first_vals))
    valid_xy <- as.data.table(xyFromCell(r[[1]], valid_cells))
    valid_xy[, cell := valid_cells]
    pts_sf <- st_as_sf(point_cells[needs_nearest], coords = c("x", "y"), crs = crs(r))
    cells_sf <- st_as_sf(valid_xy, coords = c("x", "y"), crs = crs(r))
    nearest_idx <- st_nearest_feature(pts_sf, cells_sf)
    nearest_cells <- valid_xy$cell[nearest_idx]
    point_cells$cell[which(needs_nearest)] <- nearest_cells
    pts_5070 <- st_transform(pts_sf, 5070)
    cells_5070 <- st_transform(cells_sf[nearest_idx, ], 5070)
    point_cells$cell_distance_m[which(needs_nearest)] <- as.numeric(st_distance(pts_5070, cells_5070, by_element = TRUE))
  }
  if (show_fill_progress) {
    pb <- txtProgressBar(min = 0, max = length(date_levels), style = 3)
    on.exit(close(pb), add = TRUE)
  }
  for (i in seq_along(date_levels)) {
    d <- date_levels[i]
    rows <- which(dt$date == d & is.na(dt$value))
    lyr <- which(names(r) == as.character(d))
    cells <- point_cells$cell[match(dt$zip[rows], point_cells$zip)]
    vals <- values(r[[lyr]], mat = FALSE)[cells]
    fillable <- !is.na(vals)
    if (any(fillable)) {
      idx <- rows[fillable]
      cell_dist <- point_cells$cell_distance_m[match(dt$zip[idx], point_cells$zip)]
      fill_vals <- vals[fillable]
      if (var %in% c("tmmx", "tmmn")) fill_vals <- fill_vals - 273.15
      set(dt, i = idx, j = "value", value = fill_vals)
      set(dt, i = idx, j = "value_source", value = "nearest_raster_cell")
      set(dt, i = idx, j = "fill_distance_m", value = cell_dist)
      set(dt, i = idx, j = "fill_cell", value = cells[fillable])
    }
    if (show_fill_progress) setTxtProgressBar(pb, i)
  }
  dt
}

aggregate_var_year <- function(var, year, zcta) {
  path <- download_gridmet_file(var, year)
  log_msg("Aggregating ", var, " ", year)
  r <- rast(path)
  if (is.na(crs(r))) crs(r) <- "EPSG:4326"
  r <- crop(r, ext(-125, -66, 24, 50))
  dates <- layer_dates(r, year)
  names(r) <- as.character(dates)
  zones <- make_zone_raster(r, zcta)
  z <- crop(zones$zones, r)
  names(z) <- "zip_id"
  zs <- as.data.table(zonal(r, z, fun = "mean", na.rm = TRUE))
  dt <- merge(zones$lookup, zs, by = "zip_id", all.x = TRUE, sort = FALSE)[, zip_id := NULL][]
  dt <- melt(dt, id.vars = "zip", variable.name = "date", value.name = "value", variable.factor = FALSE)
  dt[, date := as.Date(date)]
  dt <- standardize_units(dt, var)
  fill_nearest_by_layer(dt, r, zones, var)
  setnames(dt, "value", var)
  setnames(dt, "value_source", paste0(var, "_source"))
  setnames(dt, "fill_distance_m", paste0(var, "_fill_distance_m"))
  setnames(dt, "fill_cell", paste0(var, "_fill_cell"))
  dt[]
}

finalize_year <- function(dt, year) {
  if (all(c("tmmx", "tmmn") %in% names(dt))) {
    dt[, tmean_c := (tmmx + tmmn) / 2]
    dt[, temp_range_c := tmmx - tmmn]
    setnames(dt, c("tmmx", "tmmn"), c("tmax_c", "tmin_c"))
  }
  if (all(c("rmax", "rmin") %in% names(dt))) {
    dt[, rhmean_pct := (rmax + rmin) / 2]
    setnames(dt, c("rmax", "rmin"), c("rmax_pct", "rmin_pct"))
  }
  dt[, year := year]
  setcolorder(dt, c("zip", "date", "year"))
  dt
}

write_year_qc <- function(dt, year) {
  value_cols <- setdiff(names(dt), c("zip", "date", "year", grep("_(source|fill_distance_m|fill_cell)$", names(dt), value = TRUE)))
  qc <- rbindlist(lapply(value_cols, function(col) {
    data.table(
      year = year,
      variable = col,
      rows = nrow(dt),
      zips = uniqueN(dt$zip),
      dates = uniqueN(dt$date),
      missing_values = sum(is.na(dt[[col]])),
      mean = mean(dt[[col]], na.rm = TRUE),
      p01 = quantile(dt[[col]], 0.01, na.rm = TRUE),
      p50 = quantile(dt[[col]], 0.50, na.rm = TRUE),
      p99 = quantile(dt[[col]], 0.99, na.rm = TRUE)
    )
  }))
  fwrite(qc, file.path(out_dir, sprintf("gridmet_zcta_daily_qc_%04d.csv", year)))
  qc
}

process_year <- function(year, zcta) {
  dest <- file.path(out_dir, sprintf("gridmet_zcta_daily_%04d.csv.gz", year))
  if (skip_existing_year && file.exists(dest)) {
    log_msg("Skipping existing year output ", year)
    return(data.table(year = year, output = dest, skipped = TRUE))
  }

  year_dt <- NULL
  for (var in target_vars) {
    var_dt <- aggregate_var_year(var, year, zcta)
    keep_cols <- c("zip", "date", var, paste0(var, c("_source", "_fill_distance_m", "_fill_cell")))
    var_dt <- var_dt[, ..keep_cols]
    if (is.null(year_dt)) {
      year_dt <- var_dt
    } else {
      year_dt <- merge(year_dt, var_dt, by = c("zip", "date"), all = TRUE, sort = FALSE)
    }
    gc()
  }
  year_dt <- finalize_year(year_dt, year)
  setorder(year_dt, zip, date)
  fwrite(year_dt, dest)
  qc <- write_year_qc(year_dt, year)
  log_msg("Wrote ", dest, " rows=", nrow(year_dt), " missing=", sum(qc$missing_values))
  data.table(year = year, output = normalizePath(dest, winslash = "/", mustWork = TRUE), skipped = FALSE)
}

log_msg("gridMET variables: ", paste(target_vars, collapse = ", "))
log_msg("gridMET years: ", paste(range(years), collapse = "-"))
zcta <- read_zctas()
manifest <- rbindlist(lapply(years, process_year, zcta = zcta), fill = TRUE)
fwrite(manifest, file.path(out_dir, "gridmet_zcta_daily_manifest.csv"))
log_msg("Done. Outputs written to ", normalizePath(out_dir, winslash = "/"))
