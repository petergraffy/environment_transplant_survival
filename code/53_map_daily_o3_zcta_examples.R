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

o3_dir <- file.path("data", "processed", "o3_zcta_daily")
out_dir <- file.path("output", "figures", "o3_daily_zcta_maps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seed <- as.integer(Sys.getenv("DAILY_O3_MAP_SEED", "20250610"))
n_maps <- as.integer(Sys.getenv("DAILY_O3_MAP_N", "4"))
start_year <- as.integer(Sys.getenv("DAILY_O3_MAP_START_YEAR", "2005"))
end_year <- as.integer(Sys.getenv("DAILY_O3_MAP_END_YEAR", "2024"))
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

monthly_files <- file.path(
  o3_dir,
  sprintf("o3_zcta_daily_%04d_%02d.csv.gz", rep(start_year:end_year, each = 12), rep(1:12, end_year - start_year + 1))
)
missing_files <- monthly_files[!file.exists(monthly_files)]
if (length(missing_files)) stop("Missing daily ozone files:\n", paste(missing_files, collapse = "\n"), call. = FALSE)

date_pool <- seq(
  as.Date(sprintf("%04d-01-01", start_year)),
  as.Date(sprintf("%04d-12-31", end_year)),
  by = "day"
)
sample_dates <- sort(sample(date_pool, n_maps))

log_msg("Reading 2020 ZCTA boundaries")
zcta <- st_read(zcta_path, quiet = TRUE)[, c("ZCTA5CE20", "geometry")]
names(zcta)[1] <- "zip"
conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326))
zcta <- suppressWarnings(st_crop(zcta, st_transform(conus_bbox, st_crs(zcta))))
zcta <- st_transform(zcta, 5070)
zcta$zip <- sprintf("%05s", zcta$zip)

read_day <- function(day) {
  year <- as.integer(format(day, "%Y"))
  month <- as.integer(format(day, "%m"))
  path <- file.path(o3_dir, sprintf("o3_zcta_daily_%04d_%02d.csv.gz", year, month))
  log_msg("Reading ", as.character(day), " from ", basename(path))
  dt <- read_dt(path, select = c("zip", "date", "o3_ppb", "value_source"))
  dt <- dt[as.Date(date) == day]
  dt[, `:=`(zip = sprintf("%05s", zip), date = as.Date(date))]
  dt
}

daily <- rbindlist(lapply(sample_dates, read_day), fill = TRUE)

summary <- daily[, .(
  rows = .N,
  zips = uniqueN(zip),
  missing_values = sum(is.na(o3_ppb)),
  nearest_fill_values = sum(value_source == "nearest_raster_cell", na.rm = TRUE),
  min = min(o3_ppb, na.rm = TRUE),
  p01 = quantile(o3_ppb, 0.01, na.rm = TRUE),
  median = median(o3_ppb, na.rm = TRUE),
  mean = mean(o3_ppb, na.rm = TRUE),
  p99 = quantile(o3_ppb, 0.99, na.rm = TRUE),
  max = max(o3_ppb, na.rm = TRUE)
), by = date]
fwrite(summary, file.path(out_dir, "o3_zcta_daily_random_day_value_summary.csv"))

limits <- quantile(daily$o3_ppb, c(0.01, 0.99), na.rm = TRUE)

plot_one <- function(day) {
  dat <- daily[date == day, .(zip, value = o3_ppb)]
  map_data <- merge(zcta, dat, by = "zip", all.x = TRUE)
  missing_polygons <- sum(is.na(map_data$value))
  fwrite(
    data.table(date = day, map_polygons = nrow(map_data), missing_polygons = missing_polygons),
    file.path(out_dir, sprintf("o3_zcta_daily_map_missingness_%s.csv", format(day, "%Y%m%d")))
  )

  ggplot(map_data) +
    geom_sf(aes(fill = pmin(pmax(value, limits[1]), limits[2])), color = NA) +
    scale_fill_viridis(
      option = "plasma",
      direction = -1,
      name = "ppb",
      limits = limits,
      na.value = "grey90"
    ) +
    coord_sf(datum = NA) +
    labs(
      title = paste0("Daily MDA8 Ozone by ZCTA5, ", format(day, "%Y-%m-%d")),
      subtitle = "CONUS daily ozone aggregated to 2020 ZCTA5; color scale clipped to sampled-day 1st-99th percentiles"
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 15, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 10, color = "grey30", margin = margin(b = 8)),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      plot.margin = margin(8, 8, 8, 8)
    )
}

log_msg("Writing maps")
plots <- setNames(lapply(sample_dates, plot_one), format(sample_dates, "%Y%m%d"))
for (nm in names(plots)) {
  ggsave(
    filename = file.path(out_dir, sprintf("o3_zcta_daily_%s.png", nm)),
    plot = plots[[nm]],
    width = 10,
    height = 6,
    dpi = 220,
    bg = "white"
  )
}

panel <- wrap_plots(plots, ncol = 1) +
  plot_annotation(
    title = "Daily Ozone ZCTA5 Aggregation Random-Day QA Maps",
    subtitle = paste0("Random seed: ", seed),
    theme = theme(plot.title = element_text(face = "bold", size = 18))
  )
ggsave(
  filename = file.path(out_dir, "o3_zcta_daily_random_day_panel.png"),
  plot = panel,
  width = 10,
  height = 5.8 * length(plots),
  dpi = 220,
  bg = "white"
)

log_msg("Maps written to ", normalizePath(out_dir, winslash = "/"))
