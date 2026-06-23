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
  library(stringr)
  library(tibble)
})

sf_use_s2(FALSE)

analysis_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_analysis_dataset.csv.gz"
)
center_path <- file.path("output", "center_maps", "transplant_centers_mapped.csv")
zcta_path <- file.path("data", "cache", "cb_2020_us_zcta520_500k", "cb_2020_us_zcta520_500k.shp")
fig_dir <- file.path("output", "figures", "waitlist_geography_ab_maps")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(analysis_path)) stop("Missing analysis dataset: ", analysis_path, call. = FALSE)
if (!file.exists(center_path)) stop("Missing center map data: ", center_path, call. = FALSE)
if (!file.exists(zcta_path)) stop("Missing cached ZCTA shapefile: ", zcta_path, call. = FALSE)

clean_zip <- function(x) sprintf("%05s", as.character(x))

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

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

count_breaks <- c(-Inf, 49, 99, 249, 499, Inf)
count_labels <- c("<50", "50-99", "100-249", "250-499", "500+")
candidate_counts <- candidate_counts %>%
  mutate(
    count_bin = cut(
      waitlist_candidates,
      breaks = count_breaks,
      labels = count_labels,
      include.lowest = TRUE,
      right = TRUE
    )
  )

log_msg("Reading CONUS ZCTA boundaries")
zcta <- st_read(zcta_path, quiet = TRUE)[, c("ZCTA5CE20", "geometry")]
names(zcta)[1] <- "zip"
zcta$zip <- clean_zip(zcta$zip)
conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326))
zcta <- suppressWarnings(st_crop(zcta, st_transform(conus_bbox, st_crs(zcta))))
zcta <- st_transform(zcta, 5070)

candidate_map <- merge(zcta, candidate_counts, by = "zip", all.x = TRUE)
candidate_map$count_bin <- factor(candidate_map$count_bin, levels = count_labels)

candidate_palette <- c(
  "<50" = "#d7e8f4",
  "50-99" = "#a6bddb",
  "100-249" = "#67a9cf",
  "250-499" = "#2166ac",
  "500+" = "#762a83"
)

state_layer <- NULL
if (requireNamespace("maps", quietly = TRUE)) {
  state_layer <- ggplot2::map_data("state")
}

centers <- read_csv(center_path, show_col_types = FALSE) %>%
  filter(!is.na(lon), !is.na(lat), lon >= -125, lon <= -66, lat >= 24, lat <= 50) %>%
  mutate(
    dominant_organ = factor(dominant_organ, levels = c("Heart", "Kidney", "Liver", "Lung")),
    total_candidates = as.numeric(total_candidates)
  )

organ_palette <- c(
  Heart = "#0072B2",
  Kidney = "#D55E00",
  Liver = "#009E73",
  Lung = "#CC79A7"
)

theme_geo <- function() {
  theme_void(base_size = 11.5) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0, margin = margin(b = 5)),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.box = "vertical",
      plot.margin = margin(7, 8, 7, 8)
    )
}

p_candidates <- ggplot(candidate_map) +
  geom_sf(aes(fill = count_bin), color = NA) +
  scale_fill_manual(
    values = candidate_palette,
    na.value = "#f8f8f8",
    drop = FALSE,
    name = "Candidates",
    guide = guide_legend(
      nrow = 1,
      byrow = TRUE,
      keywidth = unit(18, "pt"),
      keyheight = unit(9, "pt")
    )
  ) +
  coord_sf(datum = NA) +
  labs(title = "A. Waitlist candidates by ZCTA") +
  theme_geo()

p_centers <- ggplot() +
  {
    if (!is.null(state_layer)) {
      geom_polygon(
        data = state_layer,
        aes(x = long, y = lat, group = group),
        fill = "#f4f4f1",
        color = "#c7c7c2",
        linewidth = 0.25
      )
    }
  } +
  geom_point(
    data = centers,
    aes(x = lon, y = lat, size = total_candidates, color = dominant_organ),
    alpha = 0.48
  ) +
  coord_quickmap(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
  scale_color_manual(values = organ_palette, drop = FALSE, name = "Dominant organ") +
  scale_size_continuous(
    range = c(1.2, 7.2),
    labels = scales::label_comma(),
    name = "Candidates"
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 3.2, alpha = 1)),
    size = guide_legend(order = 2, nrow = 1)
  ) +
  labs(title = "B. Transplant centers represented in cohort") +
  theme_geo()

panel <- p_candidates / p_centers + plot_layout(heights = c(1, 1))

png_path <- file.path(fig_dir, "waitlist_candidates_and_transplant_centers_ab_map.png")
pdf_path <- file.path(fig_dir, "waitlist_candidates_and_transplant_centers_ab_map.pdf")
ggsave(png_path, panel, width = 9.2, height = 11.2, dpi = 600, bg = "white")
ggsave(pdf_path, panel, width = 9.2, height = 11.2, device = grDevices::pdf, bg = "white")

write_csv(candidate_counts, file.path(fig_dir, "waitlist_candidates_by_zcta_counts.csv"))
write_csv(
  tibble(
    zctas_with_candidates = nrow(candidate_counts),
    total_distinct_candidates_across_zctas = sum(candidate_counts$waitlist_candidates),
    total_waitlist_episodes = sum(candidate_counts$waitlist_episodes),
    mapped_centers = nrow(centers)
  ),
  file.path(fig_dir, "waitlist_geography_ab_map_summary.csv")
)
write_csv(
  tibble(figure_path = c(png_path, pdf_path)),
  file.path(fig_dir, "waitlist_geography_ab_map_manifest.csv")
)

log_msg("Wrote waitlist geography A/B maps to ", normalizePath(fig_dir, winslash = "/"))
