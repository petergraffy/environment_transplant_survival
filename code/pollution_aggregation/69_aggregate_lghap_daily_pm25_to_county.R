#!/usr/bin/env Rscript
source(file.path("code", "r_runtime.R"))
ensure_user_library()
suppressPackageStartupMessages({
  library(data.table)
  library(ncdf4)
  library(Matrix)
  library(terra)
  library(sf)
  library(exactextractr)
})
data.table::setDTthreads(1)

county_out <- "data/processed/lghap_pm25_county_daily"
county_cache <- "data/cache/lghap_pm25_county_2020_area_weights_v1.rds"
county_shp <- "data/cache/cb_2020_us_county_500k/cb_2020_us_county_500k.shp"
county_script <- "code/pollution_aggregation/69_aggregate_lghap_daily_pm25_to_county.R"

county_log <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

county_archives <- function(years) {
  directories <- c("data/raw/lghap_pm25_daily", Sys.getenv("LGHAP_COUNTY_DOWNLOADS", file.path(Sys.getenv("USERPROFILE", path.expand("~")), "Downloads")))
  files <- unlist(lapply(directories, list.files,
    pattern = "^LGHAP\\.Global_PM25\\.D001\\.M[0-9]{6}\\.zip$", full.names = TRUE))
  x <- data.table(path = normalizePath(files, winslash = "/", mustWork = TRUE))
  x[, ym := sub(".*M([0-9]{6})\\.zip$", "\\1", basename(path))]
  x <- unique(x, by = "ym")
  x[, `:=`(year = as.integer(substr(ym, 1, 4)), month = as.integer(substr(ym, 5, 6)))]
  x <- x[year %in% years][order(year, month)]
  stopifnot(nrow(x) == length(years) * 12L)
  for (i in seq_len(nrow(x))) {
    first <- as.Date(sprintf("%04d-%02d-01", x$year[i], x$month[i]))
    dates <- seq(first, seq(first, by = "month", length.out = 2)[2] - 1, by = "day")
    expected <- paste0("LGHAP.Global_PM25.D001.A", format(dates, "%Y%m%d"), ".nc")
    listing <- unzip(x$path[i], list = TRUE)
    listing <- listing[grepl("\\.nc$", listing$Name), ]
    stopifnot(setequal(listing$Name, expected), !anyDuplicated(listing$Name), all(listing$Length > 0))
  }
  x
}

county_read_grid <- function(path) {
  nc <- nc_open(path)
  on.exit(nc_close(nc))
  stopifnot(identical(vapply(nc$var$PM25$dim, function(d) d$name, ""), c("lat", "lon")))
  lat <- nc$dim$lat$vals
  lon <- nc$dim$lon$vals
  stopifnot(all(diff(lat) > 0), all(diff(lon) > 0))
  iy <- which(lat >= 24 & lat <= 50)
  ix <- which(lon >= -125 & lon <= -66)
  list(lat = lat[iy], lon = lon[ix], start = c(min(iy), min(ix)), count = c(length(iy), length(ix)))
}

county_template <- function(grid) {
  dx <- median(diff(grid$lon)); dy <- median(diff(grid$lat))
  rast(nrows = length(grid$lat), ncols = length(grid$lon),
    xmin = min(grid$lon) - dx / 2, xmax = max(grid$lon) + dx / 2,
    ymin = min(grid$lat) - dy / 2, ymax = max(grid$lat) + dy / 2, crs = "EPSG:4326")
}

county_prepare <- function(archives) {
  if (file.exists(county_cache)) return(readRDS(county_cache))
  counties <- st_read(county_shp, quiet = TRUE)
  counties <- counties[!counties$STATEFP %in% c("02", "15", "60", "66", "69", "72", "78"), ]
  counties <- st_transform(st_make_valid(counties[order(counties$GEOID), ]), 4326)
  stopifnot(nrow(counties) == 3108L, !anyDuplicated(counties$GEOID))
  first_member <- unzip(archives$path[1], list = TRUE)$Name[1]
  extract_dir <- tempfile("county_grid_")
  dir.create(extract_dir)
  unzip(archives$path[1], files = first_member, exdir = extract_dir)
  grid <- county_read_grid(file.path(extract_dir, first_member))
  unlink(file.path(extract_dir, first_member))
  r <- county_template(grid)
  values(r) <- 1
  county_log("Building exact county/grid overlap areas once")
  parts <- exact_extract(r, counties, include_cell = TRUE, include_area = TRUE, progress = TRUE)
  weights <- rbindlist(lapply(seq_along(parts), function(i) {
    p <- parts[[i]]
    data.table(county_id = i, cell = p$cell, area = p$area * p$coverage_fraction)
  }))[area > 0]
  cells <- sort(unique(weights$cell))
  w <- sparseMatrix(i = weights$county_id, j = match(weights$cell, cells), x = weights$area,
    dims = c(nrow(counties), length(cells)))
  areas <- as.numeric(rowSums(w))
  stopifnot(all(areas > 0))
  w <- Diagonal(x = 1 / areas) %*% w
  # NetCDF is lat-by-lon, ascending latitude; raster cell IDs run from north to south.
  source_indices <- length(grid$lat) - (cells - 1L) %/% length(grid$lon) +
    ((cells - 1L) %% length(grid$lon)) * length(grid$lat)
  lookup <- data.table(county_fips = counties$GEOID, state_fips = counties$STATEFP,
    county_name = counties$NAME)
  cache <- list(grid = grid, cells = cells, source_indices = source_indices,
    weights = w, county_area_m2 = areas, lookup = lookup, counties = counties,
    boundary_sha256 = unname(tools::sha256sum(county_shp)), method = "exact_area_weighted_mean")
  saveRDS(cache, county_cache, compress = FALSE)
  fwrite(lookup, file.path(county_out, "county_lookup_2020_conus.csv"))
  cache
}

county_read_values <- function(path, cache) {
  nc <- nc_open(path)
  on.exit(nc_close(nc))
  g <- cache$grid
  stopifnot(identical(nc$dim$lat$vals[seq.int(g$start[1], length.out = g$count[1])], g$lat),
    identical(nc$dim$lon$vals[seq.int(g$start[2], length.out = g$count[2])], g$lon))
  vals <- ncvar_get(nc, "PM25", start = g$start, count = g$count)
  vals[vals <= -900 | vals >= 32767] <- NA_real_
  stopifnot(all(vals[is.finite(vals)] >= 0))
  vals
}

county_month <- function(job) {
  cache <- get("county_weights", envir = .GlobalEnv)
  year <- job$year; month <- job$month
  dest <- file.path(county_out, sprintf("lghap_pm25_county_daily_%04d_%02d.csv.gz", year, month))
  qc_dest <- sub("\\.csv\\.gz$", "_qc.csv", dest)
  first <- as.Date(sprintf("%04d-%02d-01", year, month))
  dates <- seq(first, seq(first, by = "month", length.out = 2)[2] - 1, by = "day")
  if (file.exists(dest) && file.exists(qc_dest)) {
    q <- fread(qc_dest)
    if (q$sha256 == unname(tools::sha256sum(dest)) && q$rows == 3108L * length(dates) && q$missing_values == 0L) {
      county_log(year, "-", sprintf("%02d", month), " verified checkpoint; skipping")
      return(q)
    }
  }
  extract_dir <- tempfile(sprintf("county_%04d_%02d_", year, month))
  dir.create(extract_dir)
  on.exit(unlink(list.files(extract_dir, full.names = TRUE)), add = TRUE)
  county_log("START ", year, "-", sprintf("%02d", month), " (", length(dates), " days)")
  result <- vector("list", length(dates))
  for (i in seq_along(dates)) {
    member <- paste0("LGHAP.Global_PM25.D001.A", format(dates[i], "%Y%m%d"), ".nc")
    unzip(job$path, files = member, exdir = extract_dir)
    path <- file.path(extract_dir, member)
    vals <- county_read_values(path, cache)
    v <- vals[cache$source_indices]
    valid <- is.finite(v)
    v[!valid] <- 0
    coverage <- as.numeric(cache$weights %*% as.numeric(valid))
    means <- as.numeric(cache$weights %*% v) / coverage
    value_source <- rep("exact_area_weighted_mean", length(means))
    fill_distance_m <- rep(NA_real_, length(means))
    fill_cell <- rep(NA_integer_, length(means))
    # Independently verify the indexed sparse calculation against raster extraction.
    if (i == 1L || any(!is.finite(means))) {
      r <- county_template(cache$grid)
      values(r) <- as.vector(t(vals[rev(seq_len(nrow(vals))), , drop = FALSE]))
      names(r) <- "pm25_ug_m3"
    }
    if (i == 1L) {
      check <- match(c("06037", "17031", "36061"), cache$lookup$county_fips)
      direct <- exact_extract(r, cache$counties[check, ], "weighted_mean", weights = "area", progress = FALSE)
      stopifnot(isTRUE(all.equal(means[check], direct, tolerance = 1e-6, check.attributes = FALSE)))
    }
    missing <- which(!is.finite(means))
    if (length(missing)) {
      # Extract all county points to avoid terra's single-point search-radius drop behavior.
      pts <- suppressWarnings(st_point_on_surface(cache$counties))
      nearest <- extract(r, vect(pts), ID = FALSE, search_radius = 500000)
      means[missing] <- nearest$pm25_ug_m3[missing]
      fill_distance_m[missing] <- nearest$distance[missing]
      fill_cell[missing] <- nearest$cell[missing]
      value_source[missing] <- "nearest_raster_cell"
      stopifnot(all(is.finite(fill_distance_m[missing])), all(fill_distance_m[missing] <= 500000))
      county_log("FILL ", dates[i], " counties=", paste(cache$lookup$county_fips[missing], collapse = ","),
        " distances_m=", paste(round(fill_distance_m[missing]), collapse = ","))
    }
    if (any(!is.finite(means))) stop("No valid county coverage on ", dates[i], ": ",
      paste(cache$lookup$county_fips[!is.finite(means)], collapse = ", "))
    stopifnot(all(means >= 0), all(coverage >= 0), all(coverage <= 1 + 1e-8))
    result[[i]] <- cbind(copy(cache$lookup), data.table(date = dates[i], year = year, month = month,
      pollutant = "pm25_lghap_daily", temporal_resolution = "daily", pm25_ug_m3 = means,
      valid_area_fraction = pmin(coverage, 1), value_source = value_source,
      fill_distance_m = fill_distance_m, fill_cell = fill_cell))
    unlink(path)
    if (i %% 10L == 0L || i == length(dates)) county_log(year, "-", sprintf("%02d", month), " day ", i, "/", length(dates))
  }
  dat <- rbindlist(result)
  setorder(dat, county_fips, date)
  stopifnot(nrow(dat) == 3108L * length(dates), !anyDuplicated(dat[, .(county_fips, date)]))
  staging <- sub("\\.csv\\.gz$", ".partial.csv.gz", dest)
  fwrite(dat, staging)
  if (!file.rename(staging, dest)) stop("Could not finalize ", dest)
  qc <- data.table(year = year, month = month, rows = nrow(dat), counties = uniqueN(dat$county_fips),
    dates = uniqueN(dat$date), missing_values = sum(!is.finite(dat$pm25_ug_m3)),
    nearest_fill_values = sum(dat$value_source == "nearest_raster_cell"),
    min_valid_area_fraction = min(dat$valid_area_fraction), mean_pm25 = mean(dat$pm25_ug_m3),
    min_pm25 = min(dat$pm25_ug_m3), max_pm25 = max(dat$pm25_ug_m3),
    file = basename(dest), sha256 = unname(tools::sha256sum(dest)))
  fwrite(qc, qc_dest)
  county_log("DONE ", year, "-", sprintf("%02d", month), " rows=", nrow(dat), " missing=0")
  qc
}

county_summarize <- function(years) {
  annual <- list(); qc <- list(); fills <- list()
  for (y in years) {
    pieces <- lapply(1:12, function(m) {
      path <- file.path(county_out, sprintf("lghap_pm25_county_daily_%04d_%02d.csv.gz", y, m))
      con <- gzfile(path, "rt")
      on.exit(close(con))
      dat <- fread(text = paste(readLines(con), collapse = "\n"), colClasses = c(county_fips = "character", state_fips = "character"))
      if (!"fill_distance_m" %in% names(dat)) {
        dat[, `:=`(fill_distance_m = NA_real_, fill_cell = NA_integer_)]
        close(con)
        on.exit(NULL)
        fwrite(dat, path)
        qpath <- sub("\\.csv\\.gz$", "_qc.csv", path)
        qrow <- fread(qpath)
        qrow[, `:=`(nearest_fill_values = sum(dat$value_source == "nearest_raster_cell"),
          sha256 = unname(tools::sha256sum(path)))]
        fwrite(qrow, qpath)
      }
      dat
    })
    d <- rbindlist(pieces, fill = TRUE)
    expected <- seq(as.Date(sprintf("%04d-01-01", y)), as.Date(sprintf("%04d-12-31", y)), by = "day")
    stopifnot(nrow(d) == 3108L * length(expected), !anyDuplicated(d[, .(county_fips, date)]),
      setequal(d$date, expected), all(is.finite(d$pm25_ug_m3)), all(grepl("^[0-9]{5}$", d$county_fips)))
    counts <- d[, .N, by = county_fips]
    stopifnot(nrow(counts) == 3108L, all(counts$N == length(expected)))
    fills[[as.character(y)]] <- d[value_source == "nearest_raster_cell"]
    annual[[as.character(y)]] <- d[, .(pm25_ug_m3 = mean(pm25_ug_m3), days = .N,
      min_valid_area_fraction = min(valid_area_fraction),
      nearest_fill_days = sum(value_source == "nearest_raster_cell")), by = .(county_fips, state_fips, county_name, year)]
    qc[[as.character(y)]] <- rbindlist(lapply(1:12, function(m)
      fread(file.path(county_out, sprintf("lghap_pm25_county_daily_%04d_%02d_qc.csv", y, m)))) , fill = TRUE)
    county_log("VALIDATED year ", y, ": ", nrow(d), " rows")
  }
  a <- rbindlist(annual); setorder(a, county_fips, year)
  q <- rbindlist(qc, fill = TRUE); setorder(q, year, month)
  q[is.na(nearest_fill_values), nearest_fill_values := 0L]
  fwrite(a, file.path(county_out, "lghap_pm25_county_annual.csv"))
  fwrite(q, file.path(county_out, "lghap_pm25_county_daily_manifest.csv"))
  fwrite(rbindlist(fills), file.path(county_out, "lghap_pm25_county_fill_audit.csv"))
  print(q[, .(counties = unique(counties), days = sum(dates), rows = sum(rows), missing_values = sum(missing_values)), by = year])
}

county_main <- function() {
  years <- seq.int(as.integer(Sys.getenv("LGHAP_COUNTY_START_YEAR", "2005")), as.integer(Sys.getenv("LGHAP_COUNTY_END_YEAR", "2024")))
  workers <- as.integer(Sys.getenv("LGHAP_COUNTY_WORKERS", "12"))
  dir.create(county_out, recursive = TRUE, showWarnings = FALSE)
  archives <- county_archives(years)
  county_log("Verified ", nrow(archives), " monthly archives; years ", min(years), "-", max(years))
  cache <- county_prepare(archives)
  stopifnot(identical(cache$boundary_sha256, unname(tools::sha256sum(county_shp))))
  if (identical(Sys.getenv("LGHAP_COUNTY_PREPARE_ONLY"), "true")) return(invisible(cache))
  cl <- parallel::makeCluster(workers, outfile = "logs/lghap_pm25_county_workers.log")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, function(root, script, cache_path) {
    setwd(root)
    source(script)
    assign("county_weights", readRDS(cache_path), envir = .GlobalEnv)
    NULL
  }, normalizePath(".", winslash = "/"), county_script, county_cache)
  for (y in years) {
    county_log("YEAR ", y, " START; ", workers, " monthly workers")
    jobs <- lapply(which(archives$year == y), function(i) as.list(archives[i]))
    parallel::parLapplyLB(cl, jobs, county_month)
    county_log("YEAR ", y, " COMPLETE")
  }
  county_summarize(years)
  county_log("ALL COMPLETE")
}

if (sys.nframe() == 0L) county_main()
