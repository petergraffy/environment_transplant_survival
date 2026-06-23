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

terraOptions(progress = as.integer(Sys.getenv("MONTHLY_NO2_TERRA_PROGRESS", "1")), memfrac = 0.75)
sf_use_s2(FALSE)

zip_path <- Sys.getenv("MONTHLY_NO2_ZIP", "C:/Users/Peter Graffy/Downloads/18919769.zip")
out_dir <- file.path("data", "processed", "no2_zcta_monthly")
cache_dir <- file.path("data", "cache")
extract_dir <- file.path(tempdir(), "monthly_no2_nc4_extract")
skip_existing <- identical(tolower(Sys.getenv("MONTHLY_NO2_SKIP_EXISTING_MONTH")), "true")
keep_nc <- identical(tolower(Sys.getenv("MONTHLY_NO2_KEEP_NC")), "true")

parse_years <- function(value, default = 2019:2025) {
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

years <- parse_years(Sys.getenv("MONTHLY_NO2_YEARS"), 2019:2025)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

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

make_zone_raster <- function(template, zcta, grid_name) {
  cache_path <- file.path(cache_dir, paste0("zcta_zone_raster_", grid_name, "_v1.tif"))
  lookup_path <- file.path(cache_dir, paste0("zcta_lookup_", grid_name, "_v1.csv"))
  points_path <- file.path(cache_dir, paste0("zcta_points_", grid_name, "_v1.csv"))
  if (file.exists(cache_path) && file.exists(lookup_path) && file.exists(points_path)) {
    return(list(
      zones = rast(cache_path),
      lookup = fread(lookup_path, colClasses = c(zip = "character")),
      points = fread(points_path, colClasses = c(zip = "character"))
    ))
  }

  log_msg("Building ZCTA zone raster for ", grid_name)
  template <- template[[1]]
  if (is.na(crs(template))) crs(template) <- "EPSG:4326"
  z_all <- st_transform(zcta, crs(template))
  conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326))
  z_all <- st_crop(z_all, st_transform(conus_bbox, st_crs(z_all)))
  z_raster <- st_crop(z_all, st_as_sfc(st_bbox(template)))
  z_vect <- vect(z_raster)
  zones <- rasterize(z_vect, template, field = "zip_id", touches = TRUE)
  pts <- st_point_on_surface(z_all)
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
  missing <- which(is.na(dt$no2_ppbv))
  if (!length(missing)) return(dt)

  point_dt <- zones$points[match(dt$zip[missing], zones$points$zip)]
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
  dt$no2_ppbv[idx] <- nearest[[value_col]][fillable]
  dt$value_source[missing] <- "nearest_raster_cell"
  dt$value_source[missing[!fillable]] <- "missing"
  dt$fill_distance_m[idx] <- nearest$distance[fillable]
  dt$fill_cell[idx] <- nearest$cell[fillable]
  dt
}

zonal_mean <- function(r, zcta) {
  if (is.na(crs(r))) crs(r) <- "EPSG:4326"
  r <- crop(r, ext(-125, -66, 24, 50))
  if (!"Surface_NO2" %in% names(r)) {
    names(r) <- "Surface_NO2"
  }
  r <- r[["Surface_NO2"]]
  r[r <= -900] <- NA
  zones <- make_zone_raster(r, zcta, "monthly_no2_conus_001deg")
  z <- crop(zones$zones, r)
  names(z) <- "zip_id"
  zs <- as.data.table(zonal(r, z, fun = "mean", na.rm = TRUE))
  setnames(zs, "Surface_NO2", "no2_ppbv")
  dt <- merge(zones$lookup, zs, by = "zip_id", all.x = TRUE, sort = FALSE)[, zip_id := NULL][]
  dt[, value_source := fifelse(is.na(no2_ppbv), "missing", "zonal_mean")]
  dt[, fill_distance_m := NA_real_]
  dt[, fill_cell := NA_integer_]
  fill_nearest_raster(dt, r, zones)
}

archive_members <- function(path) {
  listing <- unzip(path, list = TRUE)
  dt <- as.data.table(listing)
  dt <- dt[grepl("^monthly_mean_tropomi_lur_conus_surface_no2_[0-9]{6}.*\\.nc4$", Name)]
  dt[, ym := sub(".*_([0-9]{6})\\.v.*$", "\\1", Name)]
  dt[, `:=`(
    month = as.integer(substr(ym, 1, 2)),
    year = as.integer(substr(ym, 3, 6))
  )]
  dt[year %in% years][order(year, month)]
}

extract_member <- function(member) {
  unlink(file.path(extract_dir, "*"), recursive = TRUE, force = TRUE)
  utils::unzip(zip_path, files = member, exdir = extract_dir)
  file.path(extract_dir, member)
}

write_result <- function(dt, year, month) {
  dt[, `:=`(
    pollutant = "no2",
    year = as.integer(year),
    month = as.integer(month),
    temporal_resolution = "monthly"
  )]
  setcolorder(dt, c("zip", "pollutant", "year", "month", "temporal_resolution"))
  out <- file.path(out_dir, sprintf("no2_zcta_monthly_%04d_%02d.csv.gz", year, month))
  fwrite(dt, out)
  out
}

if (!file.exists(zip_path)) stop("Monthly NO2 zip not found: ", zip_path, call. = FALSE)

log_msg("Monthly NO2 archive: ", zip_path)
log_msg("Years: ", paste(years, collapse = ", "))

members <- archive_members(zip_path)
if (!nrow(members)) stop("No monthly NO2 members found in archive.", call. = FALSE)
log_msg("Monthly NO2 files in range: ", nrow(members))

zcta <- read_zctas()
manifest <- list()

for (i in seq_len(nrow(members))) {
  f <- members[i]
  out_path <- file.path(out_dir, sprintf("no2_zcta_monthly_%04d_%02d.csv.gz", f$year, f$month))
  if (skip_existing && file.exists(out_path)) {
    log_msg("Skipping existing ", f$year, "-", sprintf("%02d", f$month))
    manifest[[length(manifest) + 1L]] <- data.table(
      year = f$year, month = f$month, file = basename(out_path), path = normalizePath(out_path, winslash = "/", mustWork = TRUE)
    )
    next
  }

  log_msg("Processing NO2 ", f$year, "-", sprintf("%02d", f$month), " (", i, "/", nrow(members), ")")
  nc_path <- extract_member(f$Name)
  dt <- zonal_mean(rast(nc_path), zcta)
  out <- write_result(dt, f$year, f$month)

  qc <- dt[, .(
    rows = .N,
    zips = uniqueN(zip),
    missing_values = sum(is.na(no2_ppbv)),
    zonal_mean_values = sum(value_source == "zonal_mean", na.rm = TRUE),
    nearest_fill_values = sum(value_source == "nearest_raster_cell", na.rm = TRUE),
    mean = mean(no2_ppbv, na.rm = TRUE),
    p01 = as.numeric(quantile(no2_ppbv, 0.01, na.rm = TRUE)),
    p50 = as.numeric(quantile(no2_ppbv, 0.50, na.rm = TRUE)),
    p99 = as.numeric(quantile(no2_ppbv, 0.99, na.rm = TRUE))
  )]
  qc[, `:=`(year = f$year, month = f$month, output = normalizePath(out, winslash = "/", mustWork = TRUE))]
  setcolorder(qc, c("year", "month"))
  fwrite(qc, file.path(out_dir, sprintf("no2_zcta_monthly_qc_%04d_%02d.csv", f$year, f$month)))

  manifest[[length(manifest) + 1L]] <- data.table(
    year = f$year,
    month = f$month,
    file = basename(out),
    path = normalizePath(out, winslash = "/", mustWork = TRUE),
    rows = nrow(dt),
    zips = uniqueN(dt$zip),
    missing_values = sum(is.na(dt$no2_ppbv))
  )

  if (!keep_nc) unlink(nc_path, force = TRUE)
  rm(dt)
  gc()
}

manifest_dt <- rbindlist(manifest, fill = TRUE)
setorder(manifest_dt, year, month)
fwrite(manifest_dt, file.path(out_dir, "no2_zcta_monthly_manifest.csv"))

qc_files <- list.files(out_dir, pattern = "^no2_zcta_monthly_qc_[0-9]{4}_[0-9]{2}\\.csv$", full.names = TRUE)
if (length(qc_files)) {
  qc_all <- rbindlist(lapply(qc_files, fread), fill = TRUE)
  setorder(qc_all, year, month)
  fwrite(qc_all, file.path(out_dir, "no2_zcta_monthly_qc_all_months.csv"))
}

log_msg("Done. Outputs written to ", normalizePath(out_dir, winslash = "/"))
