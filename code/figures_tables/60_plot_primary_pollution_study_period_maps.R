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
  library(tibble)
  library(viridis)
})

sf_use_s2(FALSE)

pollution_dir <- file.path("data", "release", "air_pollution_zcta_parquet")
zcta_path <- file.path("data", "cache", "cb_2020_us_zcta520_500k", "cb_2020_us_zcta520_500k.shp")
fig_dir <- file.path("output", "figures", "primary_pollution_study_period_maps")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(zcta_path)) stop("Missing cached ZCTA shapefile: ", zcta_path, call. = FALSE)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

clean_zip <- function(x) sprintf("%05s", as.character(x))

log_msg("Reading CONUS ZCTA boundaries")
zcta <- st_read(zcta_path, quiet = TRUE)[, c("ZCTA5CE20", "geometry")]
names(zcta)[1] <- "zip"
zcta$zip <- clean_zip(zcta$zip)
conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326))
zcta <- suppressWarnings(st_crop(zcta, st_transform(conus_bbox, st_crs(zcta))))
zcta <- st_transform(zcta, 5070)

log_msg("Reading primary pollution release tables")
pm25 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_pm25_monthly_2005_2023.parquet")) %>%
  transmute(zip = clean_zip(zip), year = as.integer(year), value = pm25_ug_m3) %>%
  filter(year >= 2005, year <= 2023) %>%
  group_by(zip, year) %>%
  summarise(n_months = sum(!is.na(value)), annual_mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
  filter(n_months == 12L, is.finite(annual_mean)) %>%
  group_by(zip) %>%
  summarise(value = mean(annual_mean, na.rm = TRUE), years = n(), .groups = "drop") %>%
  filter(years == 19L)

o3 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_o3_monthly_2005_2023.parquet")) %>%
  transmute(zip = clean_zip(zip), year = as.integer(year), value = o3_ppb) %>%
  filter(year >= 2005, year <= 2023) %>%
  group_by(zip, year) %>%
  summarise(n_months = sum(!is.na(value)), annual_mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
  filter(n_months == 12L, is.finite(annual_mean)) %>%
  group_by(zip) %>%
  summarise(value = mean(annual_mean, na.rm = TRUE), years = n(), .groups = "drop") %>%
  filter(years == 19L)

no2 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_no2_annual_2005_2025.parquet")) %>%
  transmute(zip = clean_zip(zip), year = as.integer(year), value = no2) %>%
  filter(year >= 2005, year <= 2025, is.finite(value)) %>%
  group_by(zip) %>%
  summarise(value = mean(value, na.rm = TRUE), years = n(), .groups = "drop") %>%
  filter(years == 21L)

map_specs <- list(
  pm25 = list(data = pm25, title = quote(A.~Mean~PM[2.5]~concentration*","~2005-2025), legend = "ug/m3", palette = "magma", direction = -1),
  o3 = list(data = o3, title = quote(B.~Mean~O[3]~concentration*","~2005-2025), legend = "ppb", palette = "viridis", direction = -1),
  no2 = list(data = no2, title = quote(C.~Mean~NO[2]~concentration*","~2005-2025), legend = "ppb", palette = "inferno", direction = -1)
)

summary <- bind_rows(lapply(names(map_specs), function(nm) {
  dat <- map_specs[[nm]]$data
  tibble(
    pollutant = nm,
    mapped_zctas = nrow(dat),
    years_required = if_else(nm == "no2", 21L, 19L),
    min = min(dat$value, na.rm = TRUE),
    p01 = as.numeric(quantile(dat$value, 0.01, na.rm = TRUE)),
    median = median(dat$value, na.rm = TRUE),
    mean = mean(dat$value, na.rm = TRUE),
    p99 = as.numeric(quantile(dat$value, 0.99, na.rm = TRUE)),
    max = max(dat$value, na.rm = TRUE)
  )
}))
write_csv(summary, file.path(fig_dir, "primary_pollution_study_period_map_summary.csv"))

theme_map <- function() {
  theme_void(base_size = 11.5) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0, margin = margin(b = 4)),
      legend.position = "bottom",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      legend.key.width = unit(68, "pt"),
      legend.key.height = unit(10, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(6, 8, 6, 8)
    )
}

plot_map <- function(spec) {
  map_data <- merge(zcta, spec$data %>% select(zip, value), by = "zip", all.x = TRUE)
  limits <- as.numeric(quantile(map_data$value, c(0.01, 0.99), na.rm = TRUE))
  legend_breaks <- pretty(limits, n = 5)
  min_gap <- diff(limits) * 0.08
  legend_breaks <- legend_breaks[
    legend_breaks > limits[1] + min_gap &
      legend_breaks < limits[2] - min_gap
  ]
  legend_breaks <- sort(unique(c(limits, legend_breaks)))
  legend_labels <- sprintf("%.1f", legend_breaks)
  ggplot(map_data) +
    geom_sf(aes(fill = pmin(pmax(value, limits[1]), limits[2])), color = NA) +
    scale_fill_viridis_c(
      option = spec$palette,
      direction = spec$direction,
      name = spec$legend,
      limits = limits,
      breaks = legend_breaks,
      labels = legend_labels,
      na.value = "grey92",
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        barheight = unit(8, "pt"),
        barwidth = unit(170, "pt"),
        ticks.linewidth = 0.25
      )
    ) +
    coord_sf(datum = NA) +
    labs(title = spec$title) +
    theme_map()
}

plots <- lapply(map_specs, plot_map)

panel <- wrap_plots(plots, ncol = 1)

png_path <- file.path(fig_dir, "figure_pollution_study_period_maps_abc.png")
pdf_path <- file.path(fig_dir, "figure_pollution_study_period_maps_abc.pdf")

ggsave(png_path, panel, width = 8.5, height = 12.2, dpi = 600, bg = "white")
ggsave(pdf_path, panel, width = 8.5, height = 12.2, device = grDevices::pdf, bg = "white")

write_csv(
  tibble(figure_path = c(png_path, pdf_path)),
  file.path(fig_dir, "primary_pollution_study_period_map_manifest.csv")
)

log_msg("Wrote primary pollution study-period maps to ", normalizePath(fig_dir, winslash = "/"))
