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
  library(stringr)
})

in_dir <- file.path("output", "waitlist_environment_exposure")
fig_dir <- file.path("output", "figures", "waitlist_environment_exposure")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

theme_waitlist <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )
}

organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")

plot_forest <- function(dat, estimate_col, x_label, title, path) {
  plot_dat <- dat %>%
    filter(model == "combined_environment") %>%
    mutate(
      organ = recode(organ, !!!organ_labels),
      exposure = str_replace(exposure, " pct", "%"),
      exposure = factor(exposure, levels = rev(unique(exposure)))
    )

  p <- ggplot(plot_dat, aes(x = .data[[estimate_col]], y = exposure)) +
    geom_vline(xintercept = 1, linewidth = 0.35, color = "grey45") +
    geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.18, linewidth = 0.45) +
    geom_point(size = 1.8) +
    scale_x_log10() +
    facet_wrap(~organ, ncol = 2) +
    labs(title = title, x = x_label, y = NULL) +
    theme_waitlist()

  ggsave(path, p, width = 11, height = 7, dpi = 220, bg = "white")
  path
}

cox <- read_csv(file.path(in_dir, "waitlist_environment_adverse_cox_results.csv"), show_col_types = FALSE)
logistic <- read_csv(file.path(in_dir, "waitlist_environment_adverse_vs_favorable_logistic_results.csv"), show_col_types = FALSE)

paths <- c(
  plot_forest(
    cox,
    "hazard_ratio",
    "Cause-specific hazard ratio for death/deterioration delisting",
    "Waitlist-Period Environment and Adverse Waitlist Outcome",
    file.path(fig_dir, "combined_environment_adverse_cox_forest.png")
  ),
  plot_forest(
    logistic,
    "odds_ratio",
    "Odds ratio for adverse exit vs transplant/improvement",
    "Waitlist-Period Environment and Adverse vs Favorable Exit",
    file.path(fig_dir, "combined_environment_adverse_vs_favorable_forest.png")
  )
)

write_csv(tibble(figure_path = paths), file.path(fig_dir, "waitlist_environment_figure_manifest.csv"))
message("Wrote waitlist environment figures to ", fig_dir)
