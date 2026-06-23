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
  library(scales)
  library(sf)
  library(tibble)
})

sf_use_s2(FALSE)

pollution_dir <- file.path("data", "release", "air_pollution_zcta_parquet")
zcta_path <- file.path("data", "cache", "cb_2020_us_zcta520_500k", "cb_2020_us_zcta520_500k.shp")
cox_path <- file.path("output", "primary_waitlist_period_pollution_cox", "primary_waitlist_period_pollution_cox_results.csv")
fig_dir <- file.path("output", "figures", "pollution_attributable_relative_hazard_maps")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(zcta_path)) stop("Missing cached ZCTA shapefile: ", zcta_path, call. = FALSE)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

clean_zip <- function(x) sprintf("%05s", as.character(x))

pollutant_specs <- list(
  pm25 = list(
    label = "PM2.5",
    unit = 5,
    unit_label = "per 5 \u00b5g/m\u00b3",
    years = "2006-2023",
    parquet = "air_pollution_zcta_pm25_monthly_2005_2023.parquet",
    value_col = "pm25_ug_m3",
    table = "monthly",
    start_year = 2006L,
    end_year = 2023L
  ),
  o3 = list(
    label = "O3",
    unit = 10,
    unit_label = "per 10 ppb",
    years = "2006-2023",
    parquet = "air_pollution_zcta_o3_monthly_2005_2023.parquet",
    value_col = "o3_ppb",
    table = "monthly",
    start_year = 2006L,
    end_year = 2023L
  ),
  no2 = list(
    label = "NO2",
    unit = 10,
    unit_label = "per 10 ppb",
    years = "2006-2025",
    parquet = "air_pollution_zcta_no2_annual_2005_2025.parquet",
    value_col = "no2",
    table = "annual",
    start_year = 2006L,
    end_year = 2025L
  )
)
pollutant_panel_labels <- c(
  pm25 = "PM[2.5]",
  o3 = "O[3]",
  no2 = "NO[2]"
)

read_pollutant_mean <- function(name, spec) {
  dat <- read_parquet(file.path(pollution_dir, spec$parquet)) %>%
    transmute(
      zip = clean_zip(zip),
      year = as.integer(year),
      value = .data[[spec$value_col]]
    ) %>%
    filter(year >= spec$start_year, year <= spec$end_year)

  if (identical(spec$table, "monthly")) {
    dat %>%
      group_by(zip, year) %>%
      summarise(n_periods = sum(!is.na(value)), annual_mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
      filter(n_periods == 12L, is.finite(annual_mean)) %>%
      group_by(zip) %>%
      summarise(value = mean(annual_mean, na.rm = TRUE), years = n(), .groups = "drop") %>%
      filter(years == (spec$end_year - spec$start_year + 1L)) %>%
      mutate(exposure = name)
  } else {
    dat %>%
      filter(is.finite(value)) %>%
      group_by(zip) %>%
      summarise(value = mean(value, na.rm = TRUE), years = n(), .groups = "drop") %>%
      filter(years == (spec$end_year - spec$start_year + 1L)) %>%
      mutate(exposure = name)
  }
}

log_msg("Reading CONUS ZCTA boundaries")
zcta <- st_read(zcta_path, quiet = TRUE)[, c("ZCTA5CE20", "geometry")]
names(zcta)[1] <- "zip"
zcta$zip <- clean_zip(zcta$zip)
conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326))
zcta <- suppressWarnings(st_crop(zcta, st_transform(conus_bbox, st_crs(zcta))))
zcta <- st_transform(zcta, 5070)

log_msg("Reading study-period ZCTA pollution means")
pollution_means <- bind_rows(lapply(names(pollutant_specs), function(name) {
  read_pollutant_mean(name, pollutant_specs[[name]])
}))

pollution_reference <- pollution_means %>%
  group_by(exposure) %>%
  summarise(
    national_median = median(value, na.rm = TRUE),
    zctas = n(),
    .groups = "drop"
  )

log_msg("Computing model-implied relative hazards versus national median")
cox_results <- read_csv(cox_path, show_col_types = FALSE) %>%
  transmute(
    exposure,
    organ = organ_label,
    beta = log(hazard_ratio),
    primary_hr = hazard_ratio,
    primary_conf_low = conf_low,
    primary_conf_high = conf_high
  )

map_values <- pollution_means %>%
  inner_join(pollution_reference, by = "exposure") %>%
  inner_join(cox_results, by = "exposure", relationship = "many-to-many") %>%
  mutate(
    exposure_unit = vapply(exposure, function(x) pollutant_specs[[x]]$unit, numeric(1)),
    pollutant_label = vapply(exposure, function(x) pollutant_specs[[x]]$label, character(1)),
    years = vapply(exposure, function(x) pollutant_specs[[x]]$years, character(1)),
    relative_hazard = exp(beta * ((value - national_median) / exposure_unit)),
    pollutant_panel = factor(exposure, levels = names(pollutant_specs)),
    organ = factor(organ, levels = c("Heart", "Kidney", "Liver", "Lung"))
  )

map_summary <- map_values %>%
  group_by(exposure, pollutant_label, years, organ) %>%
  summarise(
    zctas = n(),
    national_median_exposure = first(national_median),
    primary_hr = first(primary_hr),
    primary_conf_low = first(primary_conf_low),
    primary_conf_high = first(primary_conf_high),
    min_relative_hazard = min(relative_hazard, na.rm = TRUE),
    p01_relative_hazard = as.numeric(quantile(relative_hazard, 0.01, na.rm = TRUE)),
    median_relative_hazard = median(relative_hazard, na.rm = TRUE),
    p99_relative_hazard = as.numeric(quantile(relative_hazard, 0.99, na.rm = TRUE)),
    max_relative_hazard = max(relative_hazard, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(map_summary, file.path(fig_dir, "pollution_attributable_relative_hazard_map_summary.csv"))

plot_dat <- merge(
  zcta,
  map_values %>% select(zip, exposure, pollutant_panel, organ, relative_hazard),
  by = "zip",
  all.x = FALSE
)

panel_limits <- plot_dat %>%
  st_drop_geometry() %>%
  group_by(exposure, organ) %>%
  summarise(
    lo = as.numeric(quantile(relative_hazard, 0.01, na.rm = TRUE)),
    hi = as.numeric(quantile(relative_hazard, 0.99, na.rm = TRUE)),
    .groups = "drop"
  )

global_limits <- panel_limits %>%
  summarise(
    lo = min(lo, na.rm = TRUE),
    hi = max(hi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    lo = floor(lo * 10) / 10,
    hi = ceiling(hi * 10) / 10
  )

legend_breaks <- seq(global_limits$lo, global_limits$hi, by = 0.1)

plot_dat <- plot_dat %>%
  mutate(
    relative_hazard_winsor = pmin(pmax(relative_hazard, global_limits$lo), global_limits$hi)
  )

write_csv(
  tibble(
    scale = "shared_panel_p01_to_p99_rounded",
    lo = global_limits$lo,
    hi = global_limits$hi,
    breaks = paste(sprintf("%.1f", legend_breaks), collapse = ", ")
  ),
  file.path(fig_dir, "pollution_attributable_relative_hazard_shared_scale.csv")
)

theme_map <- function() {
  theme_void(base_size = 11.5) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0, margin = margin(b = 4)),
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "grey94", color = "grey70", linewidth = 0.35),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10.5),
      legend.key.width = unit(130, "pt"),
      legend.key.height = unit(12, "pt"),
      panel.spacing = unit(4, "pt"),
      plot.margin = margin(5, 8, 5, 8)
    )
}

panel <- ggplot(plot_dat) +
  geom_sf(aes(fill = relative_hazard_winsor), color = NA) +
  facet_grid(
    pollutant_panel ~ organ,
    labeller = labeller(pollutant_panel = as_labeller(pollutant_panel_labels, label_parsed))
  ) +
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "#E6E6E6",
    high = "#D7191C",
    midpoint = 1,
    limits = c(global_limits$lo, global_limits$hi),
    breaks = legend_breaks,
    labels = function(x) sprintf("%.1f", x),
    oob = scales::squish,
    name = "Relative hazard",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(13, "pt"),
      barwidth = unit(520, "pt"),
      ticks = TRUE,
      ticks.linewidth = 0.3
    )
  ) +
  coord_sf(datum = NA) +
  theme_map()

png_path <- file.path(fig_dir, "figure_pollution_attributable_relative_hazard_maps.png")
pdf_path <- file.path(fig_dir, "figure_pollution_attributable_relative_hazard_maps.pdf")

ggsave(png_path, panel, width = 13.5, height = 8.2, dpi = 500, bg = "white")
ggsave(pdf_path, panel, width = 13.5, height = 8.2, device = grDevices::pdf, bg = "white")

write_csv(
  tibble(
    figure_path = c(png_path, pdf_path),
    description = "Model-implied pollutant-attributable relative hazard by ZCTA versus national median pollutant concentration."
  ),
  file.path(fig_dir, "pollution_attributable_relative_hazard_map_manifest.csv")
)

log_msg("Wrote pollution-attributable relative hazard maps to ", normalizePath(fig_dir, winslash = "/"))
