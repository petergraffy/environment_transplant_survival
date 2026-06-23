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
  library(readr)
  library(tibble)
})

curve_path <- file.path("output", "formal_waitlist_environment", "cif_by_exposure_quartile_curves.csv")
out_dir <- file.path("output", "figures", "formal_waitlist_environment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_overlay <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

curves <- read_csv(curve_path, show_col_types = FALSE) %>%
  filter(time_days <= 3650, exposure != "rmax") %>%
  mutate(
    quartile = factor(quartile, levels = paste0("Q", 1:4)),
    years = time_days / 365.25
  )

overlay <- bind_rows(
  curves %>%
    filter(outcome == "Death/deterioration delisting") %>%
    transmute(
      organ,
      organ_label,
      exposure,
      exposure_label,
      quartile,
      years,
      estimate,
      curve = "Death/deterioration CIF"
    ),
  curves %>%
    filter(outcome == "Transplant/improvement") %>%
    transmute(
      organ,
      organ_label,
      exposure,
      exposure_label,
      quartile,
      years,
      estimate = 1 - estimate,
      curve = "1 - transplant/improvement CIF"
    )
) %>%
  mutate(curve = factor(curve, levels = c("Death/deterioration CIF", "1 - transplant/improvement CIF")))

plot_path <- file.path(out_dir, "competing_cif_overlay_by_exposure_quartile.png")
p <- ggplot(overlay, aes(x = years, y = estimate, color = quartile, linetype = curve)) +
  geom_step(linewidth = 0.65, alpha = 0.95) +
  facet_grid(exposure_label ~ organ_label) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_color_brewer(palette = "Dark2", na.translate = FALSE) +
  scale_linetype_manual(values = c("solid", "22")) +
  labs(
    title = "Adverse CIF and Transplant/Improvement-Free Probability by Exposure Quartile",
    subtitle = "Solid: death/deterioration cumulative incidence. Dotted: 1 - transplant/improvement cumulative incidence. Follow-up truncated at 10 years.",
    x = "Years from waitlist start",
    y = "Probability",
    color = "Quartile",
    linetype = NULL
  ) +
  theme_overlay()

ggsave(plot_path, p, width = 15, height = 11, dpi = 240, bg = "white")

write_csv(
  tibble(
    figure = "competing_cif_overlay_by_exposure_quartile",
    path = plot_path
  ),
  file.path(out_dir, "competing_cif_overlay_manifest.csv")
)

message("Wrote competing CIF overlay to ", plot_path)
