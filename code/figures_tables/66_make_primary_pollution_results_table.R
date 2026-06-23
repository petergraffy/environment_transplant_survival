#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(dplyr)
  library(gt)
  library(readr)
  library(tidyr)
})

primary_path <- file.path("output", "primary_waitlist_period_pollution_cox", "primary_waitlist_period_pollution_cox_results.csv")
interaction_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_subgroup_cox",
  "primary_waitlist_period_pollution_subgroup_interactions.csv"
)
out_dir <- file.path("output", "tables", "primary_waitlist_period_pollution")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pollutant_levels <- c("pm25", "o3", "no2")
pollutant_labels <- c(
  pm25 = "PM\u2082.\u2085 per 5 \u00b5g/m\u00b3",
  o3 = "O\u2083 per 10 ppb",
  no2 = "NO\u2082 per 10 ppb"
)

format_ci <- function(hr, low, high) {
  sprintf("%.2f (%.2f-%.2f)", hr, low, high)
}

format_p <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

primary <- read_csv(primary_path, show_col_types = FALSE) %>%
  transmute(
    organ = organ_label,
    exposure,
    pollutant = pollutant_labels[exposure],
    `Observation window` = if_else(exposure_data_end_year == 2025, "2006-2025", "2006-2023"),
    `Candidates, n` = n,
    `Adverse events, n` = adverse_events,
    `Adjusted HR (95% CI)` = format_ci(hazard_ratio, conf_low, conf_high),
    `Primary p value` = format_p(p_value)
  )

interactions <- read_csv(interaction_path, show_col_types = FALSE) %>%
  mutate(
    interaction_col = case_when(
      subgroup == "age_group" ~ "Interaction p: age",
      subgroup == "sex" ~ "Interaction p: sex",
      subgroup == "race_group" ~ "Interaction p: race",
      TRUE ~ paste0("Interaction p: ", subgroup)
    ),
    interaction_p = format_p(interaction_p_value)
  ) %>%
  select(organ = organ_label, exposure, interaction_col, interaction_p) %>%
  pivot_wider(names_from = interaction_col, values_from = interaction_p)

table_dat <- primary %>%
  left_join(interactions, by = c("organ", "exposure")) %>%
  mutate(
    organ = factor(organ, levels = c("Heart", "Kidney", "Liver", "Lung")),
    pollutant = factor(pollutant, levels = pollutant_labels[pollutant_levels])
  ) %>%
  arrange(organ, pollutant) %>%
  select(
    Organ = organ,
    Pollutant = pollutant,
    `Observation window`,
    `Candidates, n`,
    `Adverse events, n`,
    `Adjusted HR (95% CI)`,
    `Primary p value`,
    `Interaction p: age`,
    `Interaction p: sex`,
    `Interaction p: race`
  )

csv_path <- file.path(out_dir, "table_primary_pollution_results_with_interactions.csv")
html_path <- file.path(out_dir, "table_primary_pollution_results_with_interactions.html")
caption_path <- file.path(out_dir, "table_primary_pollution_results_caption.txt")

write_csv(table_dat, csv_path)

gt_table <- table_dat %>%
  gt(groupname_col = "Organ") %>%
  tab_header(
    title = "Primary Waitlist-Period Pollution Models",
    subtitle = "Cause-specific Cox models for death or deterioration delisting"
  ) %>%
  cols_label(
    Pollutant = "Pollutant",
    `Observation window` = "Window",
    `Candidates, n` = "Candidates",
    `Adverse events, n` = "Events",
    `Adjusted HR (95% CI)` = "Adjusted HR (95% CI)",
    `Primary p value` = "P",
    `Interaction p: age` = "Age int. P",
    `Interaction p: sex` = "Sex int. P",
    `Interaction p: race` = "Race int. P"
  ) %>%
  fmt_number(columns = c(`Candidates, n`, `Adverse events, n`), decimals = 0, use_seps = TRUE) %>%
  tab_options(
    table.font.names = "Arial",
    table.font.size = px(13),
    heading.title.font.size = px(17),
    heading.subtitle.font.size = px(13),
    row_group.font.weight = "bold",
    column_labels.font.weight = "bold",
    table.width = pct(100)
  )

gtsave(gt_table, html_path)

writeLines(
  c(
    "Table caption:",
    "Primary cause-specific Cox models estimate the association between observed waitlist-period pollutant exposure and death or delisting due to deterioration.",
    "Models are fit separately by organ and pollutant and adjusted for age, sex, race, listing year, baseline organ score, and listing-center strata.",
    "Interaction p-values are likelihood-ratio tests comparing models with and without exposure-by-subgroup interaction terms."
  ),
  caption_path
)

message("Wrote primary results table to ", normalizePath(out_dir, winslash = "/"))
