#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(sf)
  library(viridis)
})

sf_use_s2(FALSE)

release_dir <- file.path("data", "release", "lghap_pm25_zcta_daily_parquet")
out_dir <- file.path("output", "figures", "lghap_pm25_daily_maps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seed <- as.integer(Sys.getenv("LGHAP_PM25_MAP_SEED", "20250605"))
n_maps <- as.integer(Sys.getenv("LGHAP_PM25_MAP_N", "3"))
set.seed(seed)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

zcta_path <- file.path("data", "cache", "cb_2020_us_zcta520_500k", "cb_2020_us_zcta520_500k.shp")
if (!file.exists(zcta_path)) stop("Missing cached ZCTA shapefile: ", zcta_path, call. = FALSE)

manifest_path <- file.path(release_dir, "lghap_pm25_zcta_daily_parquet_manifest.csv")
if (!file.exists(manifest_path)) stop("Missing release manifest: ", manifest_path, call. = FALSE)

log_msg("Reading 2020 ZCTA boundaries")
zcta <- st_read(zcta_path, quiet = TRUE)[, c("ZCTA5CE20", "geometry")]
names(zcta)[1] <- "zip"
conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326))
conus_bbox <- st_transform(conus_bbox, st_crs(zcta))
zcta <- suppressWarnings(st_crop(zcta, conus_bbox))
zcta <- st_transform(zcta, 5070)

manifest <- read_csv(manifest_path, show_col_types = FALSE)
date_pool <- bind_rows(lapply(seq_len(nrow(manifest)), function(i) {
  tibble::tibble(date = seq(as.Date(manifest$first_date[i]), as.Date(manifest$last_date[i]), by = "day"))
}))
sample_dates <- sort(sample(date_pool$date, n_maps))

read_day <- function(day) {
  year <- as.integer(format(day, "%Y"))
  parquet_path <- file.path(release_dir, sprintf("lghap_pm25_zcta_daily_%04d.parquet", year))
  if (!file.exists(parquet_path)) stop("Missing Parquet file: ", parquet_path, call. = FALSE)
  log_msg("Reading ", as.character(day), " from ", basename(parquet_path))
  read_parquet(
    parquet_path,
    col_select = c("zip", "date", "pm25_ug_m3"),
    as_data_frame = TRUE
  ) %>%
    filter(date == day) %>%
    mutate(zip = sprintf("%05s", zip))
}

daily <- bind_rows(lapply(sample_dates, read_day))

summary <- daily %>%
  group_by(date) %>%
  summarise(
    rows = n(),
    zips = n_distinct(zip),
    missing_values = sum(is.na(pm25_ug_m3)),
    min = min(pm25_ug_m3, na.rm = TRUE),
    p01 = quantile(pm25_ug_m3, 0.01, na.rm = TRUE),
    median = median(pm25_ug_m3, na.rm = TRUE),
    mean = mean(pm25_ug_m3, na.rm = TRUE),
    p99 = quantile(pm25_ug_m3, 0.99, na.rm = TRUE),
    max = max(pm25_ug_m3, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(summary, file.path(out_dir, "lghap_pm25_daily_map_value_summary.csv"))

limits <- quantile(daily$pm25_ug_m3, c(0.01, 0.99), na.rm = TRUE)

plot_one <- function(day) {
  dat <- daily %>% filter(date == day) %>% select(zip, value = pm25_ug_m3)
  map_data <- merge(zcta, dat, by = "zip", all.x = TRUE)
  missing_polygons <- sum(is.na(map_data$value))
  write_csv(
    tibble::tibble(date = day, map_polygons = nrow(map_data), missing_polygons = missing_polygons),
    file.path(out_dir, sprintf("lghap_pm25_daily_map_missingness_%s.csv", format(day, "%Y%m%d")))
  )

  ggplot(map_data) +
    geom_sf(aes(fill = pmin(pmax(value, limits[1]), limits[2])), color = NA) +
    scale_fill_viridis(
      option = "magma",
      direction = -1,
      name = "ug/m3",
      limits = limits,
      na.value = "grey90"
    ) +
    coord_sf(datum = NA) +
    labs(
      title = paste0("LGHAP Daily PM2.5 by ZCTA, ", format(day, "%Y-%m-%d")),
      subtitle = "2020 ZCTA5 aggregation; color scale clipped to sampled-day 1st-99th percentiles"
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
    filename = file.path(out_dir, sprintf("lghap_pm25_zcta_daily_%s.png", nm)),
    plot = plots[[nm]],
    width = 10,
    height = 6,
    dpi = 220,
    bg = "white"
  )
}

panel <- wrap_plots(plots, ncol = 1) +
  plot_annotation(
    title = "LGHAP Daily PM2.5 ZCTA Aggregation Random-Day QA Maps",
    subtitle = paste0("Random seed: ", seed),
    theme = theme(plot.title = element_text(face = "bold", size = 18))
  )
ggsave(
  filename = file.path(out_dir, "lghap_pm25_zcta_daily_random_day_panel.png"),
  plot = panel,
  width = 10,
  height = 5.8 * length(plots),
  dpi = 220,
  bg = "white"
)

log_msg("Maps written to ", normalizePath(out_dir, winslash = "/"))
