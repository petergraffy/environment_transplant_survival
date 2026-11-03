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
  library(tibble)
})

center_path <- file.path("output", "center_maps", "transplant_centers_by_organ_mapped.csv")
fig_dir <- file.path("output", "figures", "transplant_centers_by_organ_four_panel")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(center_path)) stop("Missing center-organ map data: ", center_path, call. = FALSE)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

organ_levels <- c("Heart", "Kidney", "Liver", "Lung")
organ_palette <- c(
  Heart = "#0072B2",
  Kidney = "#D55E00",
  Liver = "#009E73",
  Lung = "#CC79A7"
)

state_layer <- NULL
if (requireNamespace("maps", quietly = TRUE)) {
  state_layer <- ggplot2::map_data("state")
}

log_msg("Reading mapped centers by organ")
centers <- read_csv(center_path, show_col_types = FALSE) %>%
  filter(!is.na(lon), !is.na(lat), lon >= -125, lon <= -66, lat >= 24, lat <= 50) %>%
  mutate(
    organ_label = factor(organ_label, levels = organ_levels),
    panel_label = factor(
      recode(
        as.character(organ_label),
        Heart = "A. Heart",
        Kidney = "B. Kidney",
        Liver = "C. Liver",
        Lung = "D. Lung"
      ),
      levels = c("A. Heart", "B. Kidney", "C. Liver", "D. Lung")
    ),
    n_candidates = as.numeric(n_candidates)
  )

summary_tbl <- centers %>%
  group_by(organ_label) %>%
  summarise(
    mapped_centers = n_distinct(center_code),
    candidates = sum(n_candidates, na.rm = TRUE),
    median_candidates_per_center = median(n_candidates, na.rm = TRUE),
    max_candidates_per_center = max(n_candidates, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(summary_tbl, file.path(fig_dir, "transplant_centers_by_organ_four_panel_summary.csv"))

theme_center_map <- function() {
  theme_void(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold", size = 15, color = "grey8", hjust = 0, margin = margin(0, 0, 3, 0)),
      strip.background = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.box = "vertical",
      panel.spacing = unit(4, "pt"),
      plot.margin = margin(4, 5, 4, 5)
    )
}

panel <- ggplot() +
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
    aes(x = lon, y = lat, size = n_candidates, color = organ_label),
    alpha = 0.42
  ) +
  facet_wrap(vars(panel_label), ncol = 2) +
  coord_quickmap(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
  scale_color_manual(values = organ_palette, guide = "none") +
  scale_size_continuous(
    range = c(1.1, 9),
    breaks = c(100, 500, 1000, 2500, 5000, 10000),
    labels = label_comma(),
    name = "Candidates",
    guide = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(alpha = 0.55, color = "grey20")
    )
  ) +
  theme_center_map()

png_path <- file.path(fig_dir, "transplant_centers_by_organ_four_panel.png")
pdf_path <- file.path(fig_dir, "transplant_centers_by_organ_four_panel.pdf")
ggsave(png_path, panel, width = 11.5, height = 8.2, dpi = 900, bg = "white")
ggsave(pdf_path, panel, width = 11.5, height = 8.2, device = grDevices::pdf, bg = "white")

write_csv(
  tibble(
    figure_path = c(png_path, pdf_path),
    description = "Four-panel map of transplant centers represented in the cohort by organ, with point size proportional to organ-specific candidate count at each center."
  ),
  file.path(fig_dir, "transplant_centers_by_organ_four_panel_manifest.csv")
)

log_msg("Wrote four-panel center map to ", normalizePath(fig_dir, winslash = "/"))
