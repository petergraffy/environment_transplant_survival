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
  library(scales)
  library(sf)
  library(stringr)
  library(tibble)
})

sf_use_s2(FALSE)

analysis_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_analysis_dataset.csv.gz"
)
zcta_path <- file.path("data", "cache", "cb_2020_us_zcta520_500k", "cb_2020_us_zcta520_500k.shp")
risk_value_path <- file.path(
  "output",
  "figures",
  "pollution_absolute_1yr_survival_maps",
  "pollution_absolute_1yr_survival_map_values.csv"
)
fig_dir <- file.path("output", "figures", "candidates_heart_lung_absolute_risk_abc")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(analysis_path)) stop("Missing analysis dataset: ", analysis_path, call. = FALSE)
if (!file.exists(zcta_path)) stop("Missing cached ZCTA shapefile: ", zcta_path, call. = FALSE)
if (!file.exists(risk_value_path)) {
  stop(
    "Missing absolute risk map values. Run code/figures_tables/92_plot_pollution_absolute_1yr_survival_maps.R first: ",
    risk_value_path,
    call. = FALSE
  )
}

clean_zip <- function(x) {
  z <- str_extract(as.character(x), "[0-9]+")
  z <- if_else(nchar(z) > 5L, substr(z, 1, 5), z)
  if_else(is.na(z), NA_character_, str_pad(z, 5, side = "left", pad = "0"))
}

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

pollutant_panel_labels <- c(pm25 = "PM[2.5]", no2 = "NO[2]")
single_hue_risk_palette <- c(
  "#1B1B3A", "#243B6B", "#1F5F8B", "#1F7A8C", "#2A9D8F",
  "#7BC8A4", "#B7E4A8", "#E9F5A1"
)
candidate_palette <- c(
  "#08306B", "#08519C", "#2171B5", "#4292C6", "#6BAED6",
  "#9ECAE1", "#C6DBEF", "#DEEBF7", "#F7FBFF"
)

log_msg("Reading CONUS ZCTA boundaries")
zcta <- st_read(zcta_path, quiet = TRUE)[, c("ZCTA5CE20", "geometry")]
names(zcta)[1] <- "zip"
zcta$zip <- clean_zip(zcta$zip)
conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326))
zcta <- suppressWarnings(st_crop(zcta, st_transform(conus_bbox, st_crs(zcta))))
zcta <- st_transform(zcta, 5070)

log_msg("Reading candidate ZCTA counts")
candidate_counts <- read_csv(
  analysis_path,
  col_select = c(PERS_ID, candidate_zip),
  col_types = cols(PERS_ID = col_character(), candidate_zip = col_character())
) %>%
  filter(!is.na(candidate_zip), nzchar(candidate_zip)) %>%
  mutate(zip = clean_zip(candidate_zip)) %>%
  group_by(zip) %>%
  summarise(
    waitlist_candidates = n_distinct(PERS_ID),
    waitlist_episodes = n(),
    .groups = "drop"
  )

write_csv(candidate_counts, file.path(fig_dir, "waitlist_candidates_by_zcta_continuous_counts.csv"))

candidate_map <- merge(zcta, candidate_counts, by = "zip", all.x = TRUE) %>%
  mutate(waitlist_candidates = coalesce(waitlist_candidates, 0))

candidate_limits <- candidate_counts %>%
  summarise(
    lo = 0,
    hi = as.numeric(quantile(waitlist_candidates, 0.99, na.rm = TRUE)),
    max_candidates = max(waitlist_candidates, na.rm = TRUE),
    .groups = "drop"
  )
candidate_breaks <- pretty(c(candidate_limits$lo, candidate_limits$hi), n = 7)

risk_values <- read_csv(risk_value_path, show_col_types = FALSE) %>%
  filter(organ %in% c("Heart", "Lung"), exposure %in% c("pm25", "no2")) %>%
  transmute(
    zip = clean_zip(zip),
    organ,
    exposure,
    pollutant_panel = factor(exposure, levels = c("pm25", "no2")),
    adverse_risk_1yr_pct = adverse_risk_1yr * 100
  )

risk_map <- merge(zcta, risk_values, by = "zip", all.x = FALSE)

theme_abc_map <- function(base_size = 14) {
  theme_void(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 6, hjust = 0, margin = margin(b = 8, l = 40)),
      strip.text = element_text(face = "bold", size = base_size + 4, color = "grey8", margin = margin(7, 7, 7, 7)),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.55),
      legend.position = "bottom",
      legend.title = element_text(size = base_size + 1),
      legend.text = element_text(size = base_size),
      legend.key.height = unit(15, "pt"),
      panel.spacing = unit(6, "pt"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

p_candidates <- ggplot(candidate_map) +
  geom_sf(aes(fill = pmin(waitlist_candidates, candidate_limits$hi)), color = NA) +
  scale_fill_gradientn(
    colours = candidate_palette,
    limits = c(candidate_limits$lo, candidate_limits$hi),
    breaks = candidate_breaks,
    labels = label_comma(),
    oob = squish,
    na.value = "#f8f8f8",
    name = "Waitlist candidates",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(15, "pt"),
      barwidth = unit(480, "pt"),
      ticks = TRUE
    )
  ) +
  coord_sf(datum = NA) +
  labs(title = "Waitlist candidates by ZCTA") +
  theme_abc_map()

make_risk_plot <- function(organ_name) {
  organ_dat <- risk_map %>%
    filter(organ == organ_name)

  organ_limits <- organ_dat %>%
    st_drop_geometry() %>%
    summarise(
      lo = floor(as.numeric(quantile(adverse_risk_1yr_pct, 0.01, na.rm = TRUE)) * 10) / 10,
      hi = ceiling(as.numeric(quantile(adverse_risk_1yr_pct, 0.99, na.rm = TRUE)) * 10) / 10,
      .groups = "drop"
    )
  organ_breaks <- pretty(c(organ_limits$lo, organ_limits$hi), n = 7)

  ggplot(organ_dat) +
    geom_sf(
      aes(fill = pmin(pmax(adverse_risk_1yr_pct, organ_limits$lo), organ_limits$hi)),
      color = NA
    ) +
    facet_wrap(
      vars(pollutant_panel),
      nrow = 1,
      labeller = labeller(pollutant_panel = as_labeller(pollutant_panel_labels, label_parsed))
    ) +
    scale_fill_gradientn(
      colours = single_hue_risk_palette,
      limits = c(organ_limits$lo, organ_limits$hi),
      breaks = organ_breaks,
      labels = function(x) sprintf("%.1f%%", x),
      oob = squish,
      name = paste0(organ_name, " model-implied 1-year risk"),
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        barheight = unit(15, "pt"),
        barwidth = unit(500, "pt"),
        ticks = TRUE
      )
    ) +
    coord_sf(datum = NA) +
    labs(title = paste(organ_name, "candidates")) +
    theme_abc_map()
}

p_heart <- make_risk_plot("Heart")
p_lung <- make_risk_plot("Lung")

panel <- p_candidates / p_heart / p_lung +
  plot_layout(heights = c(0.95, 1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 24),
    plot.tag.position = c(0.015, 0.985)
  )

png_path <- file.path(fig_dir, "candidates_heart_lung_absolute_risk_abc.png")
pdf_path <- file.path(fig_dir, "candidates_heart_lung_absolute_risk_abc.pdf")
ggsave(png_path, panel, width = 12, height = 18, dpi = 360, device = grDevices::png, bg = "white")
pdf_path <- tryCatch(
  {
    ggsave(pdf_path, panel, width = 12, height = 18, device = grDevices::pdf, bg = "white")
    pdf_path
  },
  error = function(e) {
    fallback_pdf <- file.path(fig_dir, "candidates_heart_lung_absolute_risk_abc_updated.pdf")
    ggsave(fallback_pdf, panel, width = 12, height = 18, device = grDevices::pdf, bg = "white")
    fallback_pdf
  }
)

write_csv(
  tibble(
    figure_path = c(png_path, pdf_path),
    description = "ABC map with continuous waitlist candidate counts by ZCTA and organ-specific model-implied 1-year death/deterioration risk maps for heart and lung candidates."
  ),
  file.path(fig_dir, "candidates_heart_lung_absolute_risk_abc_manifest.csv")
)

write_csv(
  tibble(
    zctas_with_candidates = nrow(candidate_counts),
    total_distinct_candidates_across_zctas = sum(candidate_counts$waitlist_candidates),
    total_waitlist_episodes = sum(candidate_counts$waitlist_episodes),
    candidate_count_p99 = candidate_limits$hi,
    candidate_count_max = candidate_limits$max_candidates
  ),
  file.path(fig_dir, "candidates_heart_lung_absolute_risk_abc_summary.csv")
)

log_msg("Wrote ABC candidate/heart/lung map to ", normalizePath(fig_dir, winslash = "/"))
