#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

saf_source <- file.path("code", "saf_paths.R")
if (!file.exists(saf_source)) saf_source <- "saf_paths.R"
source(saf_source)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(haven)
  library(patchwork)
  library(readr)
  library(scales)
  library(stringr)
  library(survival)
  library(tibble)
  library(tidyr)
})

in_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_analysis_dataset.csv.gz"
)
fig_dir <- file.path("output", "figures", "organ_subgroup_aalen_johansen_cif")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

saf_paths <- get_saf_paths(required = TRUE, release = "q1_2026")
pubsaf_dir <- saf_paths$pubsaf_dir

pollutant_specs <- tibble(
  pollutant = c("pm25", "o3", "no2"),
  pollutant_label = c("PM[2.5]", "O[3]", "NO[2]"),
  exposure_col = c("pm25_waitlist_ug_m3", "o3_waitlist_ppb", "no2_waitlist"),
  event_col = c("event_pm25", "event_o3", "event_no2"),
  followup_col = c("followup_days_pm25", "followup_days_o3", "followup_days_no2"),
  end_year = c(2023L, 2023L, 2025L),
  exposure_end_date = as.Date(c("2023-12-31", "2023-12-31", "2025-12-31"))
)

quartile_labels <- c("Q1 lowest", "Q2", "Q3", "Q4 highest")
quartile_colors <- c(
  "Q1 lowest" = "#2166AC",
  "Q2" = "#67A9CF",
  "Q3" = "#F4A582",
  "Q4 highest" = "#B2182B"
)
timepoint_years <- c(1, 3, 5, 10)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

flag_yes <- function(x) {
  y <- str_to_upper(str_trim(as.character(x)))
  as.integer(y %in% c("1", "Y", "YES", "TRUE", "T"))
}

durable_vad_brand_codes <- c(
  202L, 205L, 206L, 207L, 208L, 209L, 210L, 212L, 213L, 214L,
  223L, 224L, 233L, 236L, 239L, 240L, 312L, 313L, 314L, 315L,
  316L, 319L, 322L, 327L, 330L, 333L, 334L
)

flag_durable_vad_brand <- function(...) {
  vals <- list(...)
  Reduce(`|`, lapply(vals, function(x) suppressWarnings(as.integer(as.character(x))) %in% durable_vad_brand_codes)) %>%
    as.integer()
}

pollutant_title <- function(pollutant_value) {
  switch(
    as.character(pollutant_value),
    "PM[2.5]" = bquote(PM[2.5]),
    "O[3]" = bquote(O[3]),
    "NO[2]" = bquote(NO[2]),
    as.character(pollutant_value)
  )
}

make_pollutant_dat <- function(dat, spec) {
  base <- dat %>%
    filter(
      listing_year_int <= spec$end_year,
      .data[[spec$followup_col]] > 0,
      is.finite(.data[[spec$exposure_col]]),
      !is.na(.data[[spec$event_col]]),
      !is.na(subgroup)
    ) %>%
    transmute(
      subgroup,
      observed_end_date,
      pollutant = spec$pollutant,
      pollutant_label = spec$pollutant_label,
      exposure = .data[[spec$exposure_col]],
      followup_days = pmax(.data[[spec$followup_col]], 0.5),
      adverse = as.integer(.data[[spec$event_col]] == 1),
      competing = as.integer(
        .data[[spec$event_col]] == 0 &
          transplant_or_improvement == 1L &
          observed_end_date <= spec$exposure_end_date
      )
    ) %>%
    mutate(
      event_status = factor(
        case_when(
          adverse == 1L ~ "adverse",
          competing == 1L ~ "transplant_or_improvement",
          TRUE ~ "censor"
        ),
        levels = c("censor", "adverse", "transplant_or_improvement")
      )
    )

  cuts <- as.numeric(quantile(base$exposure, probs = seq(0, 1, 0.25), na.rm = TRUE, type = 7))
  cuts[1] <- -Inf
  cuts[length(cuts)] <- Inf

  base %>%
    mutate(
      quartile = cut(
        exposure,
        breaks = cuts,
        labels = quartile_labels,
        include.lowest = TRUE,
        right = TRUE
      ),
      quartile = factor(quartile, levels = quartile_labels)
    )
}

tidy_aj <- function(dat) {
  fit <- survfit(Surv(followup_days / 365.25, event_status) ~ quartile, data = dat)
  s <- summary(fit)
  adverse_col <- match("adverse", s$states)
  transplant_col <- match("transplant_or_improvement", s$states)
  n_risk_start <- if (is.matrix(s$n.risk)) s$n.risk[, 1] else s$n.risk
  tibble(
    time = s$time,
    cif_adverse = as.numeric(s$pstate[, adverse_col]),
    cif_transplant_or_improvement = as.numeric(s$pstate[, transplant_col]),
    n_risk = as.numeric(n_risk_start),
    strata = as.character(s$strata)
  ) %>%
    mutate(
      quartile = sub("^quartile=", "", strata),
      quartile = factor(quartile, levels = quartile_labels),
      pollutant = first(dat$pollutant),
      pollutant_label = first(dat$pollutant_label),
      subgroup = first(dat$subgroup)
    )
}

tidy_aj_times <- function(dat, times) {
  fit <- survfit(Surv(followup_days / 365.25, event_status) ~ quartile, data = dat)
  s <- summary(fit, times = times, extend = TRUE)
  adverse_col <- match("adverse", s$states)
  transplant_col <- match("transplant_or_improvement", s$states)
  n_risk_start <- if (is.matrix(s$n.risk)) s$n.risk[, 1] else s$n.risk
  tibble(
    time = s$time,
    cif_adverse = as.numeric(s$pstate[, adverse_col]),
    cif_transplant_or_improvement = as.numeric(s$pstate[, transplant_col]),
    n_risk = as.numeric(n_risk_start),
    strata = as.character(s$strata)
  ) %>%
    mutate(
      quartile = sub("^quartile=", "", strata),
      quartile = factor(quartile, levels = quartile_labels),
      pollutant = first(dat$pollutant),
      pollutant_label = first(dat$pollutant_label),
      subgroup = first(dat$subgroup)
    )
}

theme_cif <- function() {
  theme_minimal(base_size = 17) +
    theme(
      text = element_text(color = "grey10"),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey45", linewidth = 0.55),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      axis.ticks.length = unit(4, "pt"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 17),
      legend.key.width = unit(28, "pt"),
      axis.title = element_text(size = 17),
      axis.text = element_text(size = 15, color = "grey20"),
      strip.text = element_text(face = "bold", size = 17),
      strip.background = element_rect(fill = "grey94", color = "grey65", linewidth = 0.5),
      panel.spacing.x = unit(12, "pt"),
      panel.spacing.y = unit(12, "pt"),
      plot.margin = margin(4, 8, 6, 8)
    )
}

build_aj_dataset <- function(dat) {
  bind_rows(lapply(seq_len(nrow(pollutant_specs)), function(i) {
    make_pollutant_dat(dat, pollutant_specs[i, ])
  })) %>%
    mutate(pollutant_label = factor(pollutant_label, levels = pollutant_specs$pollutant_label))
}

make_outputs <- function(plot_dat, stem, title, width, height) {
  event_summary <- plot_dat %>%
    group_by(pollutant, pollutant_label, subgroup, quartile) %>%
    summarise(
      n = n(),
      adverse_events = sum(event_status == "adverse", na.rm = TRUE),
      competing_events = sum(event_status == "transplant_or_improvement", na.rm = TRUE),
      censored = sum(event_status == "censor", na.rm = TRUE),
      median_exposure = median(exposure, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv(event_summary, file.path(fig_dir, paste0(stem, "_event_summary.csv")))

  split_dat <- split(plot_dat, list(plot_dat$pollutant, plot_dat$subgroup), drop = TRUE)
  curve_dat <- bind_rows(lapply(split_dat, tidy_aj)) %>%
    mutate(
      pollutant_label = factor(pollutant_label, levels = pollutant_specs$pollutant_label),
      subgroup = factor(subgroup, levels = levels(plot_dat$subgroup))
    )
  timepoint_dat <- bind_rows(lapply(split_dat, tidy_aj_times, times = timepoint_years)) %>%
    mutate(
      pollutant_label = factor(pollutant_label, levels = pollutant_specs$pollutant_label),
      subgroup = factor(subgroup, levels = levels(plot_dat$subgroup)),
      cif_adverse_percent = 100 * cif_adverse,
      cif_transplant_or_improvement_percent = 100 * cif_transplant_or_improvement
    )

  write_csv(curve_dat, file.path(fig_dir, paste0(stem, "_curve_data.csv")))
  write_csv(timepoint_dat, file.path(fig_dir, paste0(stem, "_cif_at_1_3_5_10_years.csv")))

  plot_for_horizon <- function(horizon) {
    breaks <- switch(
      as.character(horizon),
      `1` = c(0, 0.25, 0.5, 0.75, 1),
      `3` = seq(0, 3, 0.5),
      `5` = seq(0, 5, 1),
      `10` = seq(0, 10, 2),
      pretty(c(0, horizon), n = 5)
    )

    p <- ggplot(curve_dat %>% filter(time <= horizon), aes(x = time, y = cif_adverse, color = quartile)) +
      geom_step(linewidth = 0.85) +
      facet_grid(
        pollutant_label ~ subgroup,
        labeller = labeller(pollutant_label = as_labeller(setNames(lapply(levels(curve_dat$pollutant_label), pollutant_title), levels(curve_dat$pollutant_label))))
      ) +
      scale_color_manual(values = quartile_colors) +
      scale_x_continuous(
        breaks = breaks,
        labels = label_number(accuracy = 0.1, trim = TRUE),
        limits = c(0, horizon),
        expand = expansion(mult = c(0.01, 0.01))
      ) +
      scale_y_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0.02, 0.04))) +
      labs(
        title = paste0(horizon, "-Year Cumulative Incidence by Pollution Quartile: ", title),
        x = "Years",
        y = "Cumulative incidence"
      ) +
      theme_cif() +
      theme(
        plot.title = element_text(face = "bold", size = 23, hjust = 0, margin = margin(b = 8)),
        legend.position = "bottom"
      )

    p
  }

  plot_paths <- bind_rows(lapply(timepoint_years, function(horizon) {
    p <- plot_for_horizon(horizon)
    png_path <- file.path(fig_dir, paste0(stem, "_", horizon, "yr.png"))
    pdf_path <- file.path(fig_dir, paste0(stem, "_", horizon, "yr.pdf"))
    ggsave(png_path, p, width = width, height = height, dpi = 360, device = "png", bg = "white")
    ggsave(pdf_path, p, width = width, height = height, device = grDevices::pdf, bg = "white")
    tibble(horizon_years = horizon, figure_path = c(png_path, pdf_path))
  }))

  file.copy(file.path(fig_dir, paste0(stem, "_10yr.png")), file.path(fig_dir, paste0(stem, ".png")), overwrite = TRUE)
  file.copy(file.path(fig_dir, paste0(stem, "_10yr.pdf")), file.path(fig_dir, paste0(stem, ".pdf")), overwrite = TRUE)
  plot_paths <- bind_rows(
    plot_paths,
    tibble(
      horizon_years = 10,
      figure_path = c(file.path(fig_dir, paste0(stem, ".png")), file.path(fig_dir, paste0(stem, ".pdf")))
    )
  )
  write_csv(plot_paths, file.path(fig_dir, paste0(stem, "_manifest.csv")))
}

log_msg("Reading primary analysis dataset")
analysis_dat <- read_csv(in_path, show_col_types = FALSE) %>%
  mutate(
    listing_year_int = as.integer(as.character(listing_year)),
    index_date = as.Date(index_date),
    observed_end_date = as.Date(observed_end_date)
  )

log_msg("Reading heart LVAD and kidney dialysis variables from Q1 2026 SAF")
heart_candidate <- read_sas(
  file.path(pubsaf_dir, "cand_thor.sas7bdat"),
  col_select = any_of(c("PX_ID", "WL_ORG", "CAN_VAD1", "CAN_VAD2"))
) %>%
  filter(WL_ORG == "HR") %>%
  transmute(PX_ID, durable_lvad = flag_durable_vad_brand(CAN_VAD1, CAN_VAD2))

kidney_candidate <- read_sas(
  file.path(pubsaf_dir, "cand_kipa.sas7bdat"),
  col_select = any_of(c("PX_ID", "CAN_ON_DIAL", "CAN_DIAL", "CAN_DIAL_DT"))
) %>%
  transmute(
    PX_ID,
    kidney_on_dialysis = flag_yes(coalesce(as.character(CAN_ON_DIAL), as.character(CAN_DIAL))),
    CAN_DIAL_DT = as.Date(CAN_DIAL_DT)
  )

heart_dat <- analysis_dat %>%
  filter(WL_ORG == "HR") %>%
  left_join(heart_candidate, by = "PX_ID") %>%
  mutate(
    durable_lvad = coalesce(durable_lvad, 0L),
    subgroup = factor(
      if_else(durable_lvad == 1L, "Durable LVAD", "No durable LVAD"),
      levels = c("No durable LVAD", "Durable LVAD")
    )
  )

kidney_dat <- analysis_dat %>%
  filter(WL_ORG == "KI") %>%
  left_join(kidney_candidate, by = "PX_ID") %>%
  mutate(
    dialysis_years = as.numeric(index_date - CAN_DIAL_DT) / 365.25,
    dialysis_years = if_else(is.na(dialysis_years) | dialysis_years < 0, NA_real_, dialysis_years),
    subgroup = case_when(
      kidney_on_dialysis != 1L ~ "Not on dialysis",
      is.na(dialysis_years) ~ "Dialysis, unknown duration",
      dialysis_years < 1 ~ "Dialysis <1 y",
      dialysis_years < 3 ~ "Dialysis 1 to <3 y",
      dialysis_years < 5 ~ "Dialysis 3 to <5 y",
      TRUE ~ "Dialysis >=5 y"
    ),
    subgroup = factor(
      subgroup,
      levels = c(
        "Not on dialysis",
        "Dialysis <1 y",
        "Dialysis 1 to <3 y",
        "Dialysis 3 to <5 y",
        "Dialysis >=5 y",
        "Dialysis, unknown duration"
      )
    )
  )

write_csv(
  bind_rows(
    heart_dat %>% count(subgroup, name = "n") %>% mutate(organ = "Heart"),
    kidney_dat %>% count(subgroup, name = "n") %>% mutate(organ = "Kidney")
  ) %>% select(organ, subgroup, n),
  file.path(fig_dir, "heart_lvad_kidney_dialysis_subgroup_counts.csv")
)

log_msg("Building heart durable LVAD CIF curves")
heart_plot_dat <- build_aj_dataset(heart_dat)
make_outputs(
  heart_plot_dat,
  stem = "heart_durable_lvad_aalen_johansen_cif",
  title = "Heart Candidates by Durable LVAD Status",
  width = 17,
  height = 12
)

log_msg("Building kidney dialysis vintage CIF curves")
kidney_plot_dat <- build_aj_dataset(kidney_dat)
make_outputs(
  kidney_plot_dat,
  stem = "kidney_dialysis_vintage_aalen_johansen_cif",
  title = "Kidney Candidates by Dialysis Vintage at Listing",
  width = 22,
  height = 12
)

writeLines(
  c(
    "Figure note:",
    "Curves show nonparametric Aalen-Johansen cumulative incidence estimates for death or delisting due to deterioration, stratified by pollutant-specific waitlist-period exposure quartile.",
    "Transplant and delisting due to improvement are treated as competing events.",
    "Other waitlist exits, administrative end of follow-up, and pollutant data end are treated as censoring events.",
    "Heart strata are based on durable LVAD brand codes in CAND_THOR CAN_VAD1 or CAN_VAD2.",
    "Kidney strata are based on dialysis status and dialysis duration at listing."
  ),
  file.path(fig_dir, "organ_subgroup_aalen_johansen_cif_note.txt")
)

log_msg("Wrote organ subgroup Aalen-Johansen CIF curves to ", normalizePath(fig_dir, winslash = "/"))
