#!/usr/bin/env Rscript
source("code/pollution_aggregation/69_aggregate_lghap_daily_pm25_to_county.R")
library(arrow)

no2_inputs <- function() {
  home <- Sys.getenv("USERPROFILE", path.expand("~"))
  box <- Sys.getenv("NO2_COUNTY_ANENBERG_DIR", file.path(home, "Box", "NO2", "Anenberg et al 2022 Global NO2"))
  archive <- Sys.getenv("MONTHLY_NO2_ZIP", file.path(home, "Downloads", "18919769.zip"))
  legacy <- data.table(year = 2005:2020, month = NA_integer_, source_product = "mohegh_anenberg",
    path = file.path(box, sprintf("%04d_final_1km.nc", 2005:2020)), member = NA_character_)
  members <- unzip(archive, list = TRUE)$Name
  annual <- members[grepl("^annual_mean_tropomi_lur_conus_surface_no2_202[1-5]\\.v1\\.02\\.nc$", members)]
  modern <- data.table(year = as.integer(sub(".*_([0-9]{4})\\.v.*", "\\1", annual)), month = NA_integer_,
    source_product = "nawaz_tropomi_v2", path = archive, member = annual)
  monthly <- members[grepl("^monthly_mean_tropomi_lur_conus_surface_no2_[0-9]{6}\\.v1\\.02.*\\.nc4$", members)]
  key <- sub(".*_([0-9]{6})\\.v.*", "\\1", monthly)
  monthly <- data.table(year = as.integer(substr(key, 3, 6)), month = as.integer(substr(key, 1, 2)),
    source_product = "nawaz_tropomi_v2", path = archive, member = monthly)
  expected <- format(seq(as.Date("2018-05-01"), as.Date("2025-12-01"), by = "month"), "%Y-%m")
  stopifnot(all(file.exists(legacy$path)), identical(sort(modern$year), 2021:2025),
    identical(sort(sprintf("%04d-%02d", monthly$year, monthly$month)), expected))
  rbindlist(list(legacy, modern[order(year)], monthly[order(year, month)]))
}

no2_raster <- function(path) {
  r <- rast(path)
  stopifnot(nlyr(r) == 1L, is.lonlat(r), names(r) %in% c("Band1", "surface_no2", "Surface_NO2"))
  if (identical(names(r), "surface_no2")) {
    # Read the annual NetCDF latitude axis explicitly. GDAL reverses some v2 files.
    nc <- nc_open(path)
    on.exit(nc_close(nc))
    variable <- nc$var$surface_no2
    stopifnot(identical(vapply(variable$dim, function(d) d$name, ""), c("lon", "lat")),
      all(diff(nc$dim$lon$vals) > 0),
      max(abs(xFromCol(r, seq_len(ncol(r))) - nc$dim$lon$vals)) < 1e-5)
    lat_order <- order(nc$dim$lat$vals, decreasing = TRUE)
    stopifnot(max(abs(yFromRow(r, seq_len(nrow(r))) - nc$dim$lat$vals[lat_order])) < 1e-5)
    values(r) <- as.vector(ncvar_get(nc, "surface_no2")[, lat_order, drop = FALSE])
  }
  # Pad both products to the full county envelope so the denominator includes area outside the source grid.
  region <- ext(-125, -66, 24, 50)
  r <- crop(r, region, snap = "out")
  v <- values(r, mat = FALSE)
  v[v <= -900] <- NA_real_
  stopifnot(all(is.na(v) | (is.finite(v) & v >= 0)))
  values(r) <- v
  r <- extend(r, region, snap = "out")
  names(r) <- "no2_ppbv"
  r
}

no2_weights <- function(r, product) {
  county_cache <<- file.path("data/cache", paste0("no2_county_2020_", product, "_area_weights_v1.rds"))
  grid <- list(lon = xFromCol(r, seq_len(ncol(r))), lat = rev(yFromRow(r, seq_len(nrow(r)))))
  cache <- if (file.exists(county_cache)) readRDS(county_cache) else county_build_weights(grid)
  stopifnot(identical(grid, cache$grid),
    identical(cache$boundary_sha256, unname(tools::sha256sum(county_shp))))
  cache
}

no2_counties <- function(r, cache, verify) {
  v <- values(r, mat = FALSE)[cache$cells]
  valid <- is.finite(v)
  v[!valid] <- 0
  coverage <- as.numeric(cache$weights %*% as.numeric(valid))
  means <- as.numeric(cache$weights %*% v) / coverage
  if (verify) {
    idx <- match(c("06037", "17031", "36061", "12087", "11001"), cache$lookup$county_fips)
    direct <- exact_extract(r, cache$counties[idx, ], "weighted_mean", weights = "area", progress = FALSE)
    stopifnot(isTRUE(all.equal(means[idx], direct, tolerance = 1e-7)))
  }
  missing <- which(coverage == 0)
  dt <- copy(cache$lookup)
  dt[, `:=`(no2_ppbv = means, valid_area_fraction = pmin(coverage, 1),
    value_source = "exact_area_weighted_mean", fill_distance_m = NA_real_, fill_cell = NA_integer_)]
  if (length(missing)) {
    points <- suppressWarnings(st_point_on_surface(cache$counties))
    nearest <- extract(r, vect(points), ID = FALSE, search_radius = 500000)
    found <- is.finite(nearest$no2_ppbv[missing]) & is.finite(nearest$distance[missing]) &
      nearest$distance[missing] <= 500000
    if (!all(found)) stop("No valid NO2 grid cell within 500 km for ",
      paste(dt$county_fips[missing[!found]], collapse = ", "))
    dt[missing, `:=`(no2_ppbv = nearest$no2_ppbv[missing], value_source = "nearest_raster_cell",
      fill_distance_m = nearest$distance[missing], fill_cell = as.integer(nearest$cell[missing]))]
  }
  stopifnot(all(is.finite(dt$no2_ppbv)), all(dt$no2_ppbv >= 0),
    all(coverage >= 0 & coverage <= 1 + 1e-8))
  dt
}

no2_main_script <- "code/pollution_aggregation/72_aggregate_no2_to_county.R"
county_out <- "data/release/no2_county_parquet"

no2_process <- function(inputs) {
  caches <- list()
  rows <- vector("list", nrow(inputs))
  for (i in seq_len(nrow(inputs))) {
    input <- inputs[i]
    product <- input$source_product
    county_log(input$period_id, "/113 ", input$year,
      if (!is.na(input$month)) sprintf("-%02d", input$month), " ", product)
    extracted <- !is.na(input$member)
    path <- if (extracted) unzip(input$path, files = input$member, exdir = tempdir()) else input$path
    stopifnot(length(path) == 1L, file.exists(path))
    r <- no2_raster(path)
    if (extracted) unlink(path)
    grid_id <- if (is.na(input$month)) product else "nawaz_tropomi_v2_monthly"
    cache <- caches[[grid_id]]
    if (is.null(cache)) {
      cache <- no2_weights(r, grid_id)
      caches[[grid_id]] <- cache
    }
    stopifnot(terra::compareGeom(r, county_template(cache$grid), stopOnError = FALSE))
    rows[[i]] <- no2_counties(r, cache, verify = is.na(input$month) || input$month == 12L)
    rows[[i]][, `:=`(year = input$year, month = input$month, source_product = product,
      source_file = if (extracted) input$member else basename(path))]
  }
  rbindlist(rows)
}

no2_main <- function() {
  dir.create(county_out, recursive = TRUE, showWarnings = FALSE)
  dir.create("logs", showWarnings = FALSE)
  inputs <- no2_inputs()[, period_id := .I]
  workers <- as.integer(Sys.getenv("NO2_COUNTY_WORKERS", "4"))
  stopifnot(length(workers) == 1L, !is.na(workers), workers >= 1L)
  first <- inputs[c(1L, 17L, 22L)]
  initial <- no2_process(first)
  pending <- inputs[!period_id %in% first$period_id]
  if (workers == 1L) {
    remaining <- list(no2_process(pending))
  } else {
    cl <- parallel::makeCluster(workers, outfile = "logs/no2_county_workers.log")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(root, script) {
      setwd(root)
      source(script)
      NULL
    }, normalizePath(".", winslash = "/"), no2_main_script)
    jobs <- split(pending, rep_len(seq_len(workers), nrow(pending)))
    remaining <- parallel::parLapply(cl, jobs, no2_process)
  }
  all <- rbindlist(c(list(initial), remaining))
  setorder(all, year, month, county_fips)
  lookup <- unique(all[, .(county_fips, state_fips, county_name)])
  stopifnot(nrow(lookup) == 3108L, all(grepl("^[0-9]{5}$", lookup$county_fips)),
    !anyDuplicated(all[, .(county_fips, year, month)]),
    all(all[, .N, by = .(year, month)]$N == nrow(lookup)),
    nrow(all) == 3108L * nrow(inputs))
  manifests <- lapply(c("annual", "monthly"), function(resolution) {
    dat <- all[is.na(month) == (resolution == "annual")]
    if (resolution == "annual") dat[, month := NULL]
    file <- paste0("no2_county_", resolution, ".parquet")
    dest <- file.path(county_out, file)
    write_parquet(dat, dest, compression = "zstd", compression_level = 9, use_dictionary = TRUE)
    reopened <- as.data.table(read_parquet(dest))
    stopifnot(is.character(reopened$county_fips), isTRUE(all.equal(reopened, dat, tolerance = 0)))
    data.table(file = file, rows = nrow(dat), counties = nrow(lookup), periods = nrow(dat) / nrow(lookup),
      first_year = min(dat$year), last_year = max(dat$year),
      nearest_fills = sum(dat$value_source == "nearest_raster_cell"), missing_values = sum(is.na(dat$no2_ppbv)),
      size_bytes = file.size(dest), sha256 = unname(tools::sha256sum(dest)))
  })
  fwrite(rbindlist(manifests), file.path(county_out, "no2_county_parquet_manifest.csv"))
  fwrite(lookup, file.path(county_out, "county_lookup_2020_conus.csv"))
  fwrite(all[value_source == "nearest_raster_cell"], file.path(county_out, "no2_county_fill_audit.csv"))
  fwrite(all[, .(rows = .N, min_valid_area_fraction = min(valid_area_fraction),
    mean_no2_ppbv = mean(no2_ppbv), nearest_fills = sum(value_source == "nearest_raster_cell")),
    by = .(year, month, source_product, source_file)], file.path(county_out, "no2_county_qc.csv"))
  stopifnot(file.copy("docs/no2_county.md", file.path(county_out, "README.md"), overwrite = TRUE))
  print(rbindlist(manifests))
}

if (sys.nframe() == 0L) no2_main()
