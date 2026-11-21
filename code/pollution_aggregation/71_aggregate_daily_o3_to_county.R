#!/usr/bin/env Rscript
source("code/pollution_aggregation/69_aggregate_lghap_daily_pm25_to_county.R")

county_prefix <- "o3_county"
county_value <- "o3_ppb"
county_pollutant <- "o3_daily_mda8"
county_out <- "data/processed/o3_county_daily"
county_cache <- "data/cache/o3_county_2020_area_weights_v1.rds"
county_script <- "code/pollution_aggregation/71_aggregate_daily_o3_to_county.R"

ozone_jobs <- function(years) {
  downloads <- Sys.getenv("O3_COUNTY_DOWNLOADS", file.path(Sys.getenv("USERPROFILE", path.expand("~")), "Downloads"))
  unlist(lapply(years, function(year) {
    path <- file.path(downloads, sprintf("%04d.zip", year))
    if (!file.exists(path)) stop("Missing ozone archive: ", path)
    listing <- unzip(path, list = TRUE)
    listing <- listing[grepl("^[0-9]{8}\\.nc$", basename(listing$Name)), ]
    dates <- seq(as.Date(sprintf("%d-01-01", year)), as.Date(sprintf("%d-12-31", year)), by = "day")
    expected <- paste0(format(dates, "%Y%m%d"), ".nc")
    members <- basename(listing$Name)
    if (!setequal(members, expected) || anyDuplicated(members) || any(listing$Length <= 0)) {
      stop("Incomplete or duplicate daily ozone files in ", path)
    }
    listing <- listing[match(expected, members), ]
    lapply(1:12, function(month) {
      list(year = year, month = month, path = path, members = listing$Name[as.integer(format(dates, "%m")) == month])
    })
  }), recursive = FALSE)
}

ozone_grid <- function(nc) {
  variable <- nc$var$mda8
  stopifnot(identical(vapply(variable$dim, function(d) d$name, ""), c("lon", "lat")),
    identical(variable$units, "ppb"))
  list(lat = nc$dim$lat$vals, lon = nc$dim$lon$vals, lon_first = TRUE)
}

county_read_values <- function(path, cache) {
  nc <- nc_open(path)
  on.exit(nc_close(nc))
  stopifnot(identical(ozone_grid(nc), cache$grid))
  vals <- ncvar_get(nc, "mda8")
  stopifnot(all(is.na(vals) | (is.finite(vals) & vals >= 0)))
  vals
}

ozone_prepare <- function(job) {
  if (file.exists(county_cache)) return(readRDS(county_cache))
  path <- unzip(job$path, files = job$members[1], exdir = tempdir(), junkpaths = TRUE)
  on.exit(unlink(path))
  nc <- nc_open(path)
  grid <- tryCatch(ozone_grid(nc), finally = nc_close(nc))
  stopifnot(all(abs(diff(grid$lat) - 0.01) < 1e-5), all(abs(diff(grid$lon) - 0.01) < 1e-5))
  county_build_weights(grid)
}

ozone_main <- function() {
  years <- 2005:2024
  workers <- as.integer(Sys.getenv("O3_COUNTY_WORKERS", "12"))
  dir.create(county_out, recursive = TRUE, showWarnings = FALSE)
  dir.create("logs", showWarnings = FALSE)
  jobs <- ozone_jobs(years)
  county_log("Verified ", length(jobs), " months of daily MDA8 archives")
  cache <- ozone_prepare(jobs[[1]])
  stopifnot(identical(cache$boundary_sha256, unname(tools::sha256sum(county_shp))))
  fwrite(cache$lookup, file.path(county_out, "county_lookup_2020_conus.csv"))
  cl <- parallel::makeCluster(workers, outfile = "logs/o3_county_workers.log")
  on.exit(parallel::stopCluster(cl))
  parallel::clusterCall(cl, function(root, script, cache_path) {
    setwd(root)
    source(script)
    assign("county_weights", readRDS(cache_path), envir = .GlobalEnv)
    NULL
  }, normalizePath(".", winslash = "/"), county_script, county_cache)
  for (year in years) {
    county_log("YEAR ", year, " START")
    year_jobs <- Filter(function(job) job$year == year, jobs)
    parallel::parLapplyLB(cl, year_jobs, county_month)
    county_log("YEAR ", year, " COMPLETE")
  }
  county_summarize(years)
  county_log("ALL COMPLETE")
}

if (sys.nframe() == 0L) ozone_main()
