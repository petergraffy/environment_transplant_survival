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
o3_dir <- file.path("data", "processed", "o3_zcta_daily")
out_dir <- file.path("output", "figures", "random_daily_pm25_o3_maps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seed <- as.integer(Sys.getenv("DAILY_PM25_O3_MAP_SEED", "20260618"))
n_maps <- as.integer(Sys.getenv("DAILY_PM25_O3_MAP_N", "4"))
start_date <- as.Date(Sys.getenv("DAILY_PM25_O3_START_DATE", "2005-01-01"))
end_date <- as.Date(Sys.getenv("DAILY_PM25_O3_END_DATE", "2021-12-31"))
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

date_pool <- seq(start_date, end_date, by = "day")
sample_dates <- sort(sample(date_pool, n_maps))

missing_month_files <- function(days, dir, prefix) {
  ym <- unique(format(days, "%Y_%m"))
  paths <- file.path(dir, sprintf("%s_%s.csv.gz", prefix, ym))
  paths[!file.exists(paths)]
}

missing_pm25 <- missing_month_files(sample_dates, pm25_dir, "lghap_pm25_zcta_daily")
missing_o3 <- missing_month_files(sample_dates, o3_dir, "o3_zcta_daily")
if (length(missing_pm25)) stop("Missing PM2.5 files:\n", paste(missing_pm25, collapse = "\n"), call. = FALSE)
if (length(missing_o3)) stop("Missing ozone files:\n", paste(missing_o3, collapse = "\n"), call. = FALSE)

log_msg("Reading 2020 ZCTA boundaries")
zcta <- st_read(zcta_path, quiet = TRUE)[, c("ZCTA5CE20", "geometry")]
names(zcta)[1] <- "zip"
conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326))
zcta <- suppressWarnings(st_crop(zcta, st_transform(conus_bbox, st_crs(zcta))))
zcta <- st_transform(zcta, 5070)
zcta$zip <- sprintf("%05s", zcta$zip)

read_pm25_day <- function(day) {
  path <- file.path(pm25_dir, sprintf("lghap_pm25_zcta_daily_%04d_%02d.csv.gz", as.integer(format(day, "%Y")), as.integer(format(day, "%m"))))
  log_msg("Reading PM2.5 ", as.character(day), " from ", basename(path))
  dt <- read_dt(path, select = c("zip", "date", "pm25_ug_m3", "value_source"))
  dt <- dt[as.Date(date) == day]
  dt[, `:=`(zip = sprintf("%05s", zip), date = as.Date(date))]
  dt
}

read_o3_day <- function(day) {
  path <- file.path(o3_dir, sprintf("o3_zcta_daily_%04d_%02d.csv.gz", as.integer(format(day, "%Y")), as.integer(format(day, "%m"))))
  log_msg("Reading O3 ", as.character(day), " from ", basename(path))
  dt <- read_dt(path, select = c("zip", "date", "o3_ppb", "value_source"))
  dt <- dt[as.Date(date) == day]
  dt[, `:=`(zip = sprintf("%05s", zip), date = as.Date(date))]
  dt
}

pm25 <- rbindlist(lapply(sample_dates, read_pm25_day), fill = TRUE)
o3 <- rbindlist(lapply(sample_dates, read_o3_day), fill = TRUE)

summary <- merge(
  pm25[, .(
    pm25_rows = .N,
    pm25_zips = uniqueN(zip),
    pm25_missing = sum(is.na(pm25_ug_m3)),
    pm25_nearest_fill = sum(value_source == "nearest_raster_cell", na.rm = TRUE),
    pm25_mean = mean(pm25_ug_m3, na.rm = TRUE),
    pm25_p50 = median(pm25_ug_m3, na.rm = TRUE),
    pm25_p99 = as.numeric(quantile(pm25_ug_m3, 0.99, na.rm = TRUE))
  ), by = date],
  o3[, .(
    o3_rows = .N,
    o3_zips = uniqueN(zip),
    o3_missing = sum(is.na(o3_ppb)),
    o3_nearest_fill = sum(value_source == "nearest_raster_cell", na.rm = TRUE),
    o3_mean = mean(o3_ppb, na.rm = TRUE),
    o3_p50 = median(o3_ppb, na.rm = TRUE),
    o3_p99 = as.numeric(quantile(o3_ppb, 0.99, na.rm = TRUE))
  ), by = date],
  by = "date"
)
fwrite(summary, file.path(out_dir, "random_daily_pm25_o3_value_summary.csv"))

pm25_limits <- quantile(pm25$pm25_ug_m3, c(0.01, 0.99), na.rm = TRUE)
o3_limits <- quantile(o3$o3_ppb, c(0.01, 0.99), na.rm = TRUE)

plot_map <- function(dat, value_col, limits, title, legend_title, option) {
  map_data <- merge(zcta, dat[, .(zip, value = get(value_col))], by = "zip", all.x = TRUE)
  ggplot(map_data) +
    geom_sf(aes(fill = pmin(pmax(value, limits[1]), limits[2])), color = NA) +
    scale_fill_viridis(
      option = option,
      direction = -1,
      name = legend_title,
      limits = limits,
      na.value = "grey90"
    ) +
    coord_sf(datum = NA) +
    labs(title = title) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 12, margin = margin(b = 4)),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      plot.margin = margin(6, 6, 6, 6)
    )
}

plot_day <- function(day) {
  pm_plot <- plot_map(
    pm25[date == day],
    "pm25_ug_m3",
    pm25_limits,
    paste0("PM2.5, ", format(day, "%Y-%m-%d")),
    "ug/m3",
    "magma"
  )
  o3_plot <- plot_map(
    o3[date == day],
    "o3_ppb",
    o3_limits,
    paste0("Ozone, ", format(day, "%Y-%m-%d")),
    "ppb",
    "plasma"
  )
  pm_plot + o3_plot +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(
      subtitle = "2020 ZCTA5 aggregation; color scales clipped to sampled-day 1st-99th percentiles"
    )
}

log_msg("Writing paired maps")
plots <- setNames(lapply(sample_dates, plot_day), format(sample_dates, "%Y%m%d"))
for (nm in names(plots)) {
  ggsave(
    filename = file.path(out_dir, sprintf("pm25_o3_zcta_random_day_%s.png", nm)),
    plot = plots[[nm]],
    width = 13,
    height = 5.5,
    dpi = 220,
    bg = "white"
  )
}

panel <- wrap_plots(plots, ncol = 1) +
  plot_annotation(
    title = "Random Individual Daily PM2.5 and Ozone ZCTA Maps",
    subtitle = paste0("Random seed: ", seed, "; sampled from ", start_date, " to ", end_date),
    theme = theme(plot.title = element_text(face = "bold", size = 18))
  )
ggsave(
  filename = file.path(out_dir, "pm25_o3_zcta_random_day_panel.png"),
  plot = panel,
  width = 13,
  height = 5.8 * length(plots),
  dpi = 220,
  bg = "white"
)

log_msg("Maps written to ", normalizePath(out_dir, winslash = "/"))
