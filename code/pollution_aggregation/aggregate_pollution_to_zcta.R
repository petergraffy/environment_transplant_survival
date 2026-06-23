user_lib <- file.path(
  Sys.getenv("LOCALAPPDATA"),
  "R",
  "win-library",
  paste(R.version$major, sub("\\..*$", "", R.version$minor), sep = ".")
)
if (dir.exists(user_lib)) {
  .libPaths(unique(c(user_lib, .libPaths())))
}

library(data.table)
library(sf)
library(terra)

terraOptions(progress = 1, memfrac = 0.75)
sf_use_s2(FALSE)

parse_years <- function(value, default = 2005:2024) {
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

years <- parse_years(Sys.getenv("POLLUTION_YEARS"), 2005:2024)
box_dir <- "C:/Users/Peter Graffy/Box"
out_dir <- "output/zip_pollution"
cache_dir <- "data/cache"
skip_existing <- identical(tolower(Sys.getenv("POLLUTION_SKIP_EXISTING")), "true")
target_pollutants <- strsplit(Sys.getenv("POLLUTION_TARGETS", "pm25,o3,no2"), ",")[[1]]
target_pollutants <- trimws(tolower(target_pollutants))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

zcta_url <- "https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_us_zcta520_500k.zip"
zcta_zip <- file.path(cache_dir, basename(zcta_url))
zcta_dir <- file.path(cache_dir, "cb_2020_us_zcta520_500k")

download_zctas <- function() {
  if (!file.exists(file.path(zcta_dir, "cb_2020_us_zcta520_500k.shp"))) {
    dir.create(zcta_dir, recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(zcta_zip)) {
      message("Downloading 2020 Census ZCTA boundaries...")
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

pollutant_var <- function(r, pollutant) {
  nms <- names(r)
  if (pollutant == "pm25") return(which(nms == "PM25")[1])
  if (pollutant == "o3") return(grep("^O3", nms))
  if (pollutant == "no2") {
    idx <- which(nms %in% c("Band1", "surface_no2"))
    return(idx[1])
  }
  stop("Unknown pollutant: ", pollutant)
}

make_zone_raster <- function(template, zcta, grid_name) {
  cache_path <- file.path(cache_dir, paste0("zcta_zone_raster_", grid_name, "_v2.tif"))
  lookup_path <- file.path(cache_dir, paste0("zcta_lookup_", grid_name, "_v2.csv"))
  points_path <- file.path(cache_dir, paste0("zcta_points_", grid_name, "_v2.csv"))
  if (file.exists(cache_path) && file.exists(lookup_path) && file.exists(points_path)) {
    return(list(
      zones = rast(cache_path),
      lookup = fread(lookup_path, colClasses = c(zip = "character")),
      points = fread(points_path, colClasses = c(zip = "character"))
    ))
  }

  message("Building ZCTA zone raster for ", grid_name, "...")
  template <- template[[1]]
  if (is.na(crs(template))) crs(template) <- "EPSG:4326"
  z_all <- st_transform(zcta, crs(template))
  conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326))
  conus_bbox <- st_transform(conus_bbox, st_crs(z_all))
  z_all <- st_crop(z_all, conus_bbox)
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
  value_cols <- setdiff(names(dt), c(
    "zip",
    grep("_(source|fill_distance_m|fill_cell)$", names(dt), value = TRUE)
  ))
  for (j in seq_along(value_cols)) {
    value_col <- value_cols[j]
    missing <- which(is.na(dt[[value_col]]))
    if (!length(missing)) next

    point_dt <- zones$points[match(dt$zip[missing], zones$points$zip)]
    point_vect <- vect(point_dt[, .(x, y)], geom = c("x", "y"), crs = crs(r))
    nearest <- as.data.table(extract(
      r[[j]],
      point_vect,
      ID = FALSE,
      search_radius = search_radius
    ))
    nearest_value <- names(r[[j]])[1]
    fillable <- !is.na(nearest[[nearest_value]])

    if (any(!fillable)) {
      unresolved <- which(!fillable)
      r_ext <- ext(r[[j]])
      clamped_points <- copy(point_dt[unresolved])
      clamped_points[, x := pmin(pmax(x, xmin(r_ext)), xmax(r_ext))]
      clamped_points[, y := pmin(pmax(y, ymin(r_ext)), ymax(r_ext))]
      clamped_vect <- vect(clamped_points[, .(x, y)], geom = c("x", "y"), crs = crs(r))
      clamped_nearest <- as.data.table(extract(
        r[[j]],
        clamped_vect,
        ID = FALSE,
        search_radius = search_radius
      ))
      clamped_fillable <- !is.na(clamped_nearest[[nearest_value]])
      if (any(clamped_fillable)) {
        nearest[unresolved[clamped_fillable], (nearest_value) := clamped_nearest[[nearest_value]][clamped_fillable]]
        nearest[unresolved[clamped_fillable], distance := clamped_nearest$distance[clamped_fillable]]
        nearest[unresolved[clamped_fillable], cell := clamped_nearest$cell[clamped_fillable]]
        fillable <- !is.na(nearest[[nearest_value]])
      }
    }

    idx <- missing[fillable]
    dt[[value_col]][idx] <- nearest[[nearest_value]][fillable]
    dt[[paste0(value_col, "_source")]][missing] <- "nearest_raster_cell"
    dt[[paste0(value_col, "_source")]][missing[!fillable]] <- "missing"
    dt[[paste0(value_col, "_fill_distance_m")]][idx] <- nearest$distance[fillable]
    dt[[paste0(value_col, "_fill_cell")]][idx] <- nearest$cell[fillable]
  }
  dt
}

zonal_mean <- function(r, zcta, pollutant, grid_name) {
  if (is.na(crs(r))) crs(r) <- "EPSG:4326"
  r <- crop(r, ext(-125, -66, 24, 50))
  idx <- pollutant_var(r, pollutant)
  r <- r[[idx]]
  if (pollutant == "pm25") {
    r[r <= -900] <- NA
  }
  zones <- make_zone_raster(r, zcta, grid_name)
  z <- crop(zones$zones, r)
  names(z) <- "zip_id"
  zs <- as.data.table(zonal(r, z, fun = "mean", na.rm = TRUE))
  setnames(zs, "zip_id", "zip_id")
  dt <- merge(zones$lookup, zs, by = "zip_id", all.x = TRUE, sort = FALSE)[, zip_id := NULL][]
  value_cols <- setdiff(names(dt), "zip")
  for (value_col in value_cols) {
    dt[, (paste0(value_col, "_source")) := fifelse(is.na(get(value_col)), "missing", "zonal_mean")]
    dt[, (paste0(value_col, "_fill_distance_m")) := NA_real_]
    dt[, (paste0(value_col, "_fill_cell")) := NA_integer_]
  }
  fill_nearest_raster(dt, r, zones)
}

write_result <- function(dt, pollutant, temporal, year, month = NA_integer_) {
  dt[, pollutant := pollutant]
  dt[, year := year]
  dt[, month := month]
  dt[, temporal_resolution := temporal]
  setcolorder(dt, c("zip", "pollutant", "year", "month", "temporal_resolution"))
  out <- if (is.na(month)) {
    file.path(out_dir, sprintf("%s_%04d.csv.gz", pollutant, year))
  } else {
    file.path(out_dir, sprintf("%s_%04d_%02d.csv.gz", pollutant, year, month))
  }
  fwrite(dt, out)
  out
}

pm_files <- function() {
  files <- list.files(
    file.path(box_dir, "PM2.5", "Monthly", "Monthly"),
    pattern = "\\.nc$",
    recursive = TRUE,
    full.names = TRUE
  )
  dt <- data.table(path = normalizePath(files, winslash = "/", mustWork = TRUE))
  dt[, ym := sub(".*\\.([0-9]{6})-[0-9]{6}\\.nc$", "\\1", basename(path))]
  dt[, `:=`(year = as.integer(substr(ym, 1, 4)), month = as.integer(substr(ym, 5, 6)))]
  dt[year %in% years]
}

o3_files <- function() {
  files <- list.files(file.path(box_dir, "O3"), pattern = "\\.nc$", full.names = TRUE)
  dt <- data.table(path = normalizePath(files, winslash = "/", mustWork = TRUE))
  dt[, year := as.integer(sub(".*_([0-9]{4})\\.nc$", "\\1", basename(path)))]
  dt[year %in% years]
}

no2_files <- function() {
  root <- file.path(box_dir, "NO2", "Anenberg et al 2022 Global NO2")
  files <- list.files(root, pattern = "\\.nc$", full.names = TRUE)
  extra_files <- strsplit(Sys.getenv("POLLUTION_NO2_EXTRA_FILES", ""), ";", fixed = TRUE)[[1]]
  extra_files <- extra_files[nzchar(extra_files)]
  if (length(extra_files)) files <- c(files, extra_files)
  dt <- data.table(path = normalizePath(files, winslash = "/", mustWork = TRUE))
  dt[grepl("^[0-9]{4}_final_1km\\.nc$", basename(path)), year := as.integer(substr(basename(path), 1, 4))]
  dt[grepl("annual_mean_.*_([0-9]{4})\\.v", basename(path)), year := as.integer(sub(".*_([0-9]{4})\\.v.*$", "\\1", basename(path)))]
  dt[year %in% years][order(year)]
}

process_pm <- function(zcta) {
  files <- pm_files()[order(year, month)]
  message("PM2.5 monthly files in range: ", nrow(files))
  out <- character()
  for (i in seq_len(nrow(files))) {
    f <- files[i]
    dest <- file.path(out_dir, sprintf("pm25_%04d_%02d.csv.gz", f$year, f$month))
    if (skip_existing && file.exists(dest)) next
    message("PM2.5 ", f$year, "-", sprintf("%02d", f$month))
    dt <- zonal_mean(rast(f$path), zcta, "pm25", "pm25_001deg")
    setnames(
      dt,
      c("PM25", "PM25_source", "PM25_fill_distance_m", "PM25_fill_cell"),
      c("pm25_ug_m3", "value_source", "fill_distance_m", "fill_cell")
    )
    out <- c(out, write_result(dt, "pm25", "monthly", f$year, f$month))
  }
  invisible(out)
}

process_o3 <- function(zcta) {
  files <- o3_files()[order(year)]
  message("O3 annual files with monthly layers in range: ", nrow(files))
  out <- character()
  for (i in seq_len(nrow(files))) {
    f <- files[i]
    if (skip_existing && all(file.exists(file.path(out_dir, sprintf("o3_%04d_%02d.csv.gz", f$year, 1:12))))) next
    message("O3 ", f$year)
    dt <- zonal_mean(rast(f$path), zcta, "o3", "o3_001deg")
    value_cols <- grep("^O3_[0-9]+$", names(dt), value = TRUE)
    for (m in seq_along(value_cols)) {
      value_col <- value_cols[m]
      monthly <- dt[, .(
        zip,
        o3_ppb = get(value_col),
        value_source = get(paste0(value_col, "_source")),
        fill_distance_m = get(paste0(value_col, "_fill_distance_m")),
        fill_cell = get(paste0(value_col, "_fill_cell"))
      )]
      out <- c(out, write_result(monthly, "o3", "monthly", f$year, m))
    }
  }
  invisible(out)
}

process_no2 <- function(zcta) {
  files <- no2_files()[order(year)]
  message("NO2 annual files in range: ", nrow(files))
  out <- character()
  for (i in seq_len(nrow(files))) {
    f <- files[i]
    dest <- file.path(out_dir, sprintf("no2_%04d.csv.gz", f$year))
    if (skip_existing && file.exists(dest)) next
    grid <- if (grepl("_final_1km\\.nc$", basename(f$path))) "no2_global_1km" else "no2_conus_001deg"
    message("NO2 ", f$year)
    dt <- zonal_mean(rast(f$path), zcta, "no2", grid)
    value_col <- setdiff(names(dt), grep("_(source|fill_distance_m|fill_cell)$", names(dt), value = TRUE))
    value_col <- setdiff(value_col, "zip")
    setnames(
      dt,
      c(value_col, paste0(value_col, "_source"), paste0(value_col, "_fill_distance_m"), paste0(value_col, "_fill_cell")),
      c("no2", "value_source", "fill_distance_m", "fill_cell")
    )
    out <- c(out, write_result(dt, "no2", "annual", f$year))
  }
  invisible(out)
}

combine_outputs <- function() {
  files <- list.files(out_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)
  files <- files[!grepl("all_|manifest", basename(files))]
  manifest <- data.table(path = normalizePath(files, winslash = "/", mustWork = TRUE))
  manifest[, file := basename(path)]
  count_rows <- function(path) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    max(length(readLines(con, warn = FALSE)) - 1L, 0L)
  }
  read_gz_csv <- function(path) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    as.data.table(read.csv(con, check.names = FALSE, colClasses = c(zip = "character")))
  }
  manifest[, rows := vapply(path, count_rows, integer(1))]
  fwrite(manifest, file.path(out_dir, "manifest.csv"))

  by_pollutant <- split(files, sub("_.*$", "", basename(files)))
  for (nm in names(by_pollutant)) {
    message("Combining ", nm, " outputs...")
    dt <- rbindlist(lapply(by_pollutant[[nm]], read_gz_csv), fill = TRUE)
    setorder(dt, zip, year, month)
    fwrite(dt, file.path(out_dir, paste0("all_", nm, "_zip.csv.gz")))
  }
}

zcta <- read_zctas()
if ("pm25" %in% target_pollutants) process_pm(zcta)
if ("o3" %in% target_pollutants) process_o3(zcta)
if ("no2" %in% target_pollutants) process_no2(zcta)
combine_outputs()

message("Done. Outputs written to ", normalizePath(out_dir, winslash = "/"))
