#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(sf)
  library(tidyr)
  library(viridis)
})

sf_use_s2(FALSE)

map_year <- as.integer(Sys.getenv("GRIDMET_MAP_YEAR", "2024"))
gridmet_dir <- file.path("data", "processed", "gridmet_zcta_daily")
out_dir <- file.path("output", "figures", "gridmet_zcta_maps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

variables <- c(
  "tmax_c", "tmin_c", "tmean_c", "temp_range_c",
  "rmax_pct", "rmin_pct", "rhmean_pct",
  "pr", "sph", "srad", "vs", "pet", "etr", "erc"
)

variable_meta <- tibble::tribble(
  ~variable, ~label, ~legend, ~palette,
  "tmax_c", "Max temp", "deg C", "inferno",
  "tmin_c", "Min temp", "deg C", "inferno",
  "tmean_c", "Mean temp", "deg C", "inferno",
  "temp_range_c", "Temp range", "deg C", "magma",
  "rmax_pct", "Max RH", "%", "viridis",
  "rmin_pct", "Min RH", "%", "viridis",
  "rhmean_pct", "Mean RH", "%", "viridis",
  "pr", "Precipitation", "mm/day", "mako",
  "sph", "Specific humidity", "kg/kg", "viridis",
  "srad", "Solar radiation", "gridMET units", "plasma",
  "vs", "Wind speed", "m/s", "cividis",
  "pet", "Potential ET", "mm/day", "magma",
  "etr", "Reference ET", "mm/day", "magma",
  "erc", "Energy release component", "index", "rocket"
)

zcta_path <- file.path("data", "cache", "cb_2020_us_zcta520_500k", "cb_2020_us_zcta520_500k.shp")
if (!file.exists(zcta_path)) stop("Missing cached ZCTA shapefile: ", zcta_path, call. = FALSE)

log_msg("Reading 2020 ZCTA boundaries")
zcta <- st_read(zcta_path, quiet = TRUE)[, c("ZCTA5CE20", "geometry")]
names(zcta)[1] <- "zip"
conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326))
conus_bbox <- st_transform(conus_bbox, st_crs(zcta))
zcta <- suppressWarnings(st_crop(zcta, conus_bbox))
zcta <- st_transform(zcta, 5070)

gridmet_path <- file.path(gridmet_dir, sprintf("gridmet_zcta_daily_%04d.csv.gz", map_year))
if (!file.exists(gridmet_path)) stop("Missing gridMET file: ", gridmet_path, call. = FALSE)

log_msg("Reading daily gridMET values for ", map_year)
daily <- read_csv(
  gridmet_path,
  col_select = all_of(c("zip", variables)),
  col_types = cols(zip = col_character(), .default = col_double()),
  show_col_types = FALSE,
  progress = TRUE
)

log_msg("Summarising annual daily means by ZCTA")
annual <- daily %>%
  group_by(zip) %>%
  summarise(across(all_of(variables), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
rm(daily)
gc()

annual_long <- annual %>%
  pivot_longer(-zip, names_to = "variable", values_to = "value") %>%
  left_join(variable_meta, by = "variable")

write_csv(annual_long, file.path(out_dir, sprintf("gridmet_zcta_%04d_annual_means_long.csv.gz", map_year)))

plot_one <- function(variable_name) {
  meta <- variable_meta %>% filter(variable == variable_name) %>% slice(1)
  dat <- annual %>%
    select(zip, value = all_of(variable_name))
  map_data <- merge(zcta, dat, by = "zip", all.x = TRUE)
  ggplot(map_data) +
    geom_sf(aes(fill = value), color = NA) +
    scale_fill_viridis(
      option = meta$palette,
      direction = 1,
      name = meta$legend,
      na.value = "grey92"
    ) +
    coord_sf(datum = NA) +
    labs(
      title = paste0(meta$label, ", ", map_year),
      subtitle = "Annual mean of daily gridMET values by 2020 ZCTA"
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

log_msg("Writing individual maps")
plots <- setNames(lapply(variables, plot_one), variables)
for (nm in names(plots)) {
  ggsave(
    filename = file.path(out_dir, sprintf("gridmet_%s_zcta_map_%04d.png", nm, map_year)),
    plot = plots[[nm]],
    width = 10,
    height = 6,
    dpi = 220,
    bg = "white"
  )
}

log_msg("Writing panel maps")
temp_panel <- (plots$tmax_c | plots$tmin_c) / (plots$tmean_c | plots$temp_range_c) +
  plot_annotation(title = sprintf("gridMET Temperature ZCTA Maps, %04d", map_year))
humidity_panel <- (plots$rmax_pct | plots$rmin_pct) / (plots$rhmean_pct | plots$sph) +
  plot_annotation(title = sprintf("gridMET Humidity ZCTA Maps, %04d", map_year))
other_panel <- (plots$pr | plots$srad) / (plots$vs | plots$pet) / (plots$etr | plots$erc) +
  plot_annotation(title = sprintf("gridMET Hydrometeorology and Fire Context ZCTA Maps, %04d", map_year))

ggsave(file.path(out_dir, sprintf("gridmet_temperature_panel_%04d.png", map_year)), temp_panel, width = 14, height = 9, dpi = 220, bg = "white")
ggsave(file.path(out_dir, sprintf("gridmet_humidity_panel_%04d.png", map_year)), humidity_panel, width = 14, height = 9, dpi = 220, bg = "white")
ggsave(file.path(out_dir, sprintf("gridmet_other_panel_%04d.png", map_year)), other_panel, width = 14, height = 13, dpi = 220, bg = "white")

map_check <- annual_long %>%
  group_by(variable, label) %>%
  summarise(
    mapped_zips = sum(!is.na(value)),
    missing_values = sum(is.na(value)),
    min = min(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(map_check, file.path(out_dir, sprintf("gridmet_map_value_summary_%04d.csv", map_year)))

polygon_check <- lapply(variables, function(v) {
  dat <- annual %>% select(zip, value = all_of(v))
  map_data <- merge(zcta, dat, by = "zip", all.x = TRUE)
  tibble::tibble(
    variable = v,
    map_polygons = nrow(map_data),
    missing_polygons = sum(is.na(map_data$value))
  )
}) %>%
  bind_rows()
write_csv(polygon_check, file.path(out_dir, sprintf("gridmet_map_missingness_check_%04d.csv", map_year)))

log_msg("Maps written to ", normalizePath(out_dir, winslash = "/"))
