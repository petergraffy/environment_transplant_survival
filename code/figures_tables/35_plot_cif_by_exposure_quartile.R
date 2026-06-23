#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(cmprsk)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

in_path <- file.path("output", "formal_waitlist_environment", "formal_waitlist_environment_analysis_dataset.csv.gz")
out_dir <- file.path("output", "figures", "formal_waitlist_environment")
table_dir <- file.path("output", "formal_waitlist_environment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

theme_cif <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

exposure_specs <- tribble(
  ~exposure, ~term, ~label,
  "tmax", "tmax_c_waitlist_mean_exact_daily", "Tmax",
  "rmax", "rmax_pct_waitlist_mean_exact_daily", "Maximum relative humidity",
  "pm25", "pm25_waitlist_ug_m3", "PM2.5",
  "o3", "o3_waitlist_ppb", "Ozone",
  "no2", "no2_waitlist", "NO2"
)

make_quartile <- function(x) {
  qs <- quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE, names = FALSE)
  qs <- unique(qs)
  if (length(qs) < 3L) return(factor(rep(NA_character_, length(x))))
  cut(x, breaks = qs, include.lowest = TRUE, labels = paste0("Q", seq_len(length(qs) - 1L)))
}

tidy_cuminc <- function(ci, org, exposure_name, exposure_label) {
  nm <- names(ci)
  nm <- nm[!nm %in% "Tests"]
  bind_rows(lapply(nm, function(name) {
    obj <- ci[[name]]
    parts <- str_match(name, "^(.*) ([0-9]+)$")
    group <- parts[, 2]
    cause <- parts[, 3]
    tibble(
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      exposure = exposure_name,
      exposure_label = exposure_label,
      quartile = group,
      cause_code = cause,
      outcome = recode(
        cause,
        `1` = "Death/deterioration delisting",
        `2` = "Transplant/improvement",
        `3` = "Other removal",
        .default = paste("Cause", cause)
      ),
      time_days = obj$time,
      estimate = obj$est,
      variance = obj$var
    )
  }))
}

plot_cif <- function(dat, outcome_name, title, path) {
  plot_dat <- dat %>%
    filter(outcome == outcome_name, time_days <= 3650) %>%
    mutate(
      quartile = factor(quartile, levels = paste0("Q", 1:4)),
      years = time_days / 365.25
    )

  p <- ggplot(plot_dat, aes(x = years, y = estimate, color = quartile)) +
    geom_step(linewidth = 0.7, alpha = 0.95) +
    facet_grid(exposure_label ~ organ_label, scales = "free_y") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_color_brewer(palette = "Dark2", na.translate = FALSE) +
    labs(
      title = title,
      subtitle = "Nonparametric cumulative incidence by organ-specific exposure quartile; follow-up truncated at 10 years for display",
      x = "Years from waitlist start",
      y = "Cumulative incidence",
      color = "Quartile"
    ) +
    theme_cif()

  ggsave(path, p, width = 15, height = 11, dpi = 240, bg = "white")
  path
}

log_msg("Reading formal analysis dataset")
base <- read_csv(in_path, show_col_types = FALSE) %>%
  mutate(
    fstatus = case_when(
      event_type == "adverse" ~ 1L,
      event_type == "transplant_or_improvement" ~ 2L,
      event_type == "other_exit" ~ 3L,
      TRUE ~ 0L
    )
  ) %>%
  filter(WL_ORG %in% names(organ_labels), followup_days > 0)

log_msg("Computing CIF curves by exposure quartile")
cif_curves <- bind_rows(lapply(seq_len(nrow(exposure_specs)), function(i) {
  spec <- exposure_specs[i, ]
  term <- spec$term
  bind_rows(lapply(names(organ_labels), function(org) {
    dat <- base %>%
      filter(WL_ORG == org, !is.na(.data[[term]])) %>%
      mutate(exposure_quartile = make_quartile(.data[[term]])) %>%
      filter(!is.na(exposure_quartile))
    if (nrow(dat) == 0L || n_distinct(dat$exposure_quartile) < 2L) return(tibble())
    log_msg("CIF ", org, " | ", spec$label, " n=", nrow(dat), " adverse=", sum(dat$fstatus == 1L))
    ci <- cuminc(
      ftime = dat$followup_days,
      fstatus = dat$fstatus,
      group = dat$exposure_quartile,
      cencode = 0
    )
    tidy_cuminc(ci, org, spec$exposure, spec$label)
  }))
}))

write_csv(cif_curves, file.path(table_dir, "cif_by_exposure_quartile_curves.csv"))

manifest <- tibble(
  figure = c(
    "adverse_cif_by_exposure_quartile",
    "transplant_improvement_cif_by_exposure_quartile"
  ),
  path = c(
    plot_cif(
      cif_curves,
      "Death/deterioration delisting",
      "Cumulative Incidence of Death/Deterioration Delisting by Exposure Quartile",
      file.path(out_dir, "adverse_cif_by_exposure_quartile.png")
    ),
    plot_cif(
      cif_curves,
      "Transplant/improvement",
      "Cumulative Incidence of Transplant/Improvement by Exposure Quartile",
      file.path(out_dir, "transplant_improvement_cif_by_exposure_quartile.png")
    )
  )
)

write_csv(manifest, file.path(out_dir, "cif_figure_manifest.csv"))
log_msg("Wrote CIF figures to ", out_dir)
