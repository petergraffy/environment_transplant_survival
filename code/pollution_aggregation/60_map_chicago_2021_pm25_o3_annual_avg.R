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
out_dir <- file.path("output", "figures", "chicago_2021_pm25_o3_annual")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

target_year <- as.integer(Sys.getenv("CHICAGO_PM25_O3_YEAR", "2021"))

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

pm25_files <- file.path(pm25_dir, sprintf("lghap_pm25_zcta_daily_%04d_%02d.csv.gz", target_year, 1:12))
o3_files <- file.path(o3_dir, sprintf("o3_zcta_daily_%04d_%02d.csv.gz", target_year, 1:12))
missing_files <- c(pm25_files[!file.exists(pm25_files)], o3_files[!file.exists(o3_files)])
if (length(missing_files)) stop("Missing daily exposure files:\n", paste(missing_files, collapse = "\n"), call. = FALSE)

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

read_year <- function(files, value_col) {
  rbindlist(lapply(files, function(path) {
    log_msg("Reading ", basename(path))
    dt <- read_dt(path, select = c("zip", "date", value_col, "value_source"))
    dt[, `:=`(zip = sprintf("%05s", zip), date = as.Date(date))]
    dt[zip %in% zcta_chicago$zip]
  }), fill = TRUE)
}

pm25_daily <- read_year(pm25_files, "pm25_ug_m3")
o3_daily <- read_year(o3_files, "o3_ppb")

pm25_annual <- pm25_daily[, .(
  pm25_ug_m3 = mean(pm25_ug_m3, na.rm = TRUE),
  days = uniqueN(date),
  missing_values = sum(is.na(pm25_ug_m3)),
  nearest_fill_days = sum(value_source == "nearest_raster_cell", na.rm = TRUE)
), by = zip]

o3_annual <- o3_daily[, .(
  o3_ppb = mean(o3_ppb, na.rm = TRUE),
  days = uniqueN(date),
  missing_values = sum(is.na(o3_ppb)),
  nearest_fill_days = sum(value_source == "nearest_raster_cell", na.rm = TRUE)
), by = zip]

summary <- merge(
  pm25_annual[, .(
    zip,
    pm25_ug_m3,
    pm25_days = days,
    pm25_missing_values = missing_values,
    pm25_nearest_fill_days = nearest_fill_days
  )],
  o3_annual[, .(
    zip,
    o3_ppb,
    o3_days = days,
    o3_missing_values = missing_values,
    o3_nearest_fill_days = nearest_fill_days
  )],
  by = "zip",
  all = TRUE
)
fwrite(summary, file.path(out_dir, sprintf("chicago_pm25_o3_annual_average_%04d_by_zcta.csv", target_year)))

plot_one <- function(dat, value_col, title, legend_title, option) {
  map_data <- merge(zcta_chicago, dat[, .(zip, value = get(value_col))], by = "zip", all.x = TRUE)
  missing_polygons <- sum(is.na(map_data$value))
  fwrite(
    data.table(year = target_year, pollutant = value_col, map_polygons = nrow(map_data), missing_polygons = missing_polygons),
    file.path(out_dir, sprintf("chicago_%s_annual_average_%04d_missingness.csv", value_col, target_year))
  )
  limits <- quantile(map_data$value, c(0.02, 0.98), na.rm = TRUE)

  ggplot(map_data) +
    geom_sf(aes(fill = pmin(pmax(value, limits[1]), limits[2])), color = "white", linewidth = 0.12) +
    scale_fill_viridis(
      option = option,
      direction = -1,
      name = legend_title,
      limits = limits,
      na.value = "grey90"
    ) +
    coord_sf(datum = NA) +
    labs(
      title = title,
      subtitle = sprintf("%04d daily mean across %d Chicago-area ZCTAs", target_year, nrow(map_data))
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 13, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 9, color = "grey30", margin = margin(b = 8)),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      plot.margin = margin(8, 8, 8, 8)
    )
}

pm_plot <- plot_one(pm25_annual, "pm25_ug_m3", "Annual Average PM2.5", "ug/m3", "magma")
o3_plot <- plot_one(o3_annual, "o3_ppb", "Annual Average Ozone", "ppb", "plasma")

panel <- pm_plot + o3_plot +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = sprintf("Chicago ZCTA5 Annual Average PM2.5 and Ozone, %04d", target_year),
    subtitle = "Daily ZCTA estimates averaged within ZIP over the calendar year",
    theme = theme(plot.title = element_text(face = "bold", size = 17))
  )

ggsave(
  file.path(out_dir, sprintf("chicago_pm25_o3_annual_average_%04d.png", target_year)),
  panel,
  width = 13,
  height = 6.5,
  dpi = 240,
  bg = "white"
)

log_msg("Map written to ", normalizePath(out_dir, winslash = "/"))
