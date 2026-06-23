#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(sf)
  library(viridis)
})

sf_use_s2(FALSE)

pm25_dir <- file.path("data", "processed", "lghap_pm25_zcta_daily")
out_dir <- file.path("output", "figures", "chicago_pm25_zcta_random_days")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

target_year <- as.integer(Sys.getenv("CHICAGO_PM25_MAP_YEAR", "2021"))
n_maps <- as.integer(Sys.getenv("CHICAGO_PM25_MAP_N", "4"))
seed <- as.integer(Sys.getenv("CHICAGO_PM25_MAP_SEED", "20250606"))
set.seed(seed)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

read_dt <- function(path, select = NULL) {
  if (grepl("\\.gz$", path) && !requireNamespace("R.utils", quietly = TRUE)) {
    dt <- as.data.table(read_csv(path, show_col_types = FALSE, progress = FALSE))
    if (!is.null(select)) dt <- dt[, ..select]
    return(dt)
  }
  fread(path, select = select)
}

zcta_path <- file.path("data", "cache", "cb_2020_us_zcta520_500k", "cb_2020_us_zcta520_500k.shp")
if (!file.exists(zcta_path)) stop("Missing cached ZCTA shapefile: ", zcta_path, call. = FALSE)

monthly_files <- file.path(pm25_dir, sprintf("lghap_pm25_zcta_daily_%04d_%02d.csv.gz", target_year, 1:12))
missing_files <- monthly_files[!file.exists(monthly_files)]
if (length(missing_files)) stop("Missing PM2.5 files:\n", paste(missing_files, collapse = "\n"), call. = FALSE)

date_pool <- seq(as.Date(sprintf("%04d-01-01", target_year)), as.Date(sprintf("%04d-12-31", target_year)), by = "day")
sample_dates <- sort(sample(date_pool, n_maps))

log_msg("Reading 2020 ZCTA boundaries")
zcta <- st_read(zcta_path, quiet = TRUE)[, c("ZCTA5CE20", "geometry")]
names(zcta)[1] <- "zip"

chicago_bbox <- st_as_sfc(st_bbox(
  c(xmin = -88.05, ymin = 41.45, xmax = -87.35, ymax = 42.15),
  crs = 4326
))
zcta_chicago <- suppressWarnings(st_crop(zcta, st_transform(chicago_bbox, st_crs(zcta))))
zcta_chicago <- st_transform(zcta_chicago, 26916)
zcta_chicago$zip <- sprintf("%05s", zcta_chicago$zip)

read_day <- function(day) {
  month <- as.integer(format(day, "%m"))
  path <- file.path(pm25_dir, sprintf("lghap_pm25_zcta_daily_%04d_%02d.csv.gz", target_year, month))
  log_msg("Reading ", as.character(day), " from ", basename(path))
  dt <- read_dt(path, select = c("zip", "date", "pm25_ug_m3", "value_source"))
  dt <- dt[as.Date(date) == day]
  dt[, `:=`(zip = sprintf("%05s", zip), date = as.Date(date))]
  dt
}

daily <- rbindlist(lapply(sample_dates, read_day), fill = TRUE)
daily_chicago <- daily[zip %in% zcta_chicago$zip]

summary <- daily_chicago[, .(
  rows = .N,
  zips = uniqueN(zip),
  missing_values = sum(is.na(pm25_ug_m3)),
  nearest_fill_values = sum(value_source == "nearest_raster_cell", na.rm = TRUE),
  min = min(pm25_ug_m3, na.rm = TRUE),
  p01 = quantile(pm25_ug_m3, 0.01, na.rm = TRUE),
  median = median(pm25_ug_m3, na.rm = TRUE),
  mean = mean(pm25_ug_m3, na.rm = TRUE),
  p99 = quantile(pm25_ug_m3, 0.99, na.rm = TRUE),
  max = max(pm25_ug_m3, na.rm = TRUE)
), by = date]
fwrite(summary, file.path(out_dir, "chicago_pm25_zcta_random_day_value_summary.csv"))

limits <- quantile(daily_chicago$pm25_ug_m3, c(0.02, 0.98), na.rm = TRUE)

plot_one <- function(day) {
  dat <- daily_chicago[date == day, .(zip, value = pm25_ug_m3)]
  map_data <- merge(zcta_chicago, dat, by = "zip", all.x = TRUE)
  missing_polygons <- sum(is.na(map_data$value))
  fwrite(
    data.table(date = day, map_polygons = nrow(map_data), missing_polygons = missing_polygons),
    file.path(out_dir, sprintf("chicago_pm25_zcta_random_day_missingness_%s.csv", format(day, "%Y%m%d")))
  )

  ggplot(map_data) +
    geom_sf(aes(fill = pmin(pmax(value, limits[1]), limits[2])), color = "white", linewidth = 0.12) +
    scale_fill_viridis(
      option = "magma",
      direction = -1,
      name = "ug/m3",
      limits = limits,
      na.value = "grey90"
    ) +
    coord_sf(datum = NA) +
    labs(
      title = paste0("Chicago PM2.5 by ZCTA5, ", format(day, "%Y-%m-%d")),
      subtitle = sprintf("LGHAP native daily PM2.5 aggregated to 2020 ZCTA5; %d ZCTAs shown", nrow(map_data))
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 9.5, color = "grey30", margin = margin(b = 8)),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      plot.margin = margin(8, 8, 8, 8)
    )
}

log_msg("Writing Chicago maps")
plots <- setNames(lapply(sample_dates, plot_one), format(sample_dates, "%Y%m%d"))
for (nm in names(plots)) {
  ggsave(
    filename = file.path(out_dir, sprintf("chicago_pm25_zcta_%s.png", nm)),
    plot = plots[[nm]],
    width = 7.5,
    height = 6.5,
    dpi = 240,
    bg = "white"
  )
}

panel <- wrap_plots(plots, ncol = 2) +
  plot_annotation(
    title = sprintf("Chicago ZCTA5 Daily PM2.5 Random-Day QA Maps, %04d", target_year),
    subtitle = paste0("Random seed: ", seed),
    theme = theme(plot.title = element_text(face = "bold", size = 17))
  )
ggsave(
  filename = file.path(out_dir, sprintf("chicago_pm25_zcta_random_day_panel_%04d.png", target_year)),
  plot = panel,
  width = 15,
  height = 13,
  dpi = 240,
  bg = "white"
)

log_msg("Maps written to ", normalizePath(out_dir, winslash = "/"))
