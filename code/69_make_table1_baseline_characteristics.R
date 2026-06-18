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
  library(stringr)
  library(tibble)
  library(tidyr)
})

in_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_analysis_dataset.csv.gz"
)
out_dir <- file.path("output", "tables", "table1_baseline_waitlist_characteristics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")
organ_levels <- unname(organ_labels)
sex_levels <- c("Female", "Male", "Unknown")
race_levels <- c("White", "Black", "Asian", "Other/Unknown")
outcome_levels <- c(
  "transplant_or_improvement",
  "death_or_deterioration_delist",
  "other_exit",
  "administrative_censor"
)
outcome_labels <- c(
  transplant_or_improvement = "Transplant or delisting due to improvement",
  death_or_deterioration_delist = "Death or delisting due to deterioration",
  other_exit = "Other waitlist exit",
  administrative_censor = "Administrative censoring"
)

fmt_int <- function(x) {
  formatC(as.integer(round(x)), format = "d", big.mark = ",")
}

fmt_n_pct <- function(n, denom) {
  if (is.na(n) || is.na(denom) || denom == 0) return("")
  paste0(fmt_int(n), " (", sprintf("%.1f", 100 * n / denom), "%)")
}

fmt_median_iqr <- function(x, digits = 1) {
  x <- x[is.finite(x)]
  if (!length(x)) return("")
  qs <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE, type = 7)
  paste0(
    sprintf(paste0("%.", digits, "f"), qs[2]),
    " (",
    sprintf(paste0("%.", digits, "f"), qs[1]),
    " to ",
    sprintf(paste0("%.", digits, "f"), qs[3]),
    ")"
  )
}

fmt_median_iqr_range <- function(x, digits = 1) {
  x <- x[is.finite(x)]
  if (!length(x)) return("")
  qs <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE, type = 7)
  paste0(
    sprintf(paste0("%.", digits, "f"), qs[2]),
    " (",
    sprintf(paste0("%.", digits, "f"), qs[1]),
    " to ",
    sprintf(paste0("%.", digits, "f"), qs[3]),
    "); ",
    sprintf(paste0("%.", digits, "f"), min(x, na.rm = TRUE)),
    " to ",
    sprintf(paste0("%.", digits, "f"), max(x, na.rm = TRUE))
  )
}

component_cont_row <- function(dat, org, var, characteristic, row_order, digits = 1, transform = identity) {
  dat %>%
    filter(WL_ORG == org) %>%
    mutate(.component_value = transform(.data[[var]])) %>%
    group_by(organ) %>%
    summarise(value = fmt_median_iqr(.component_value, digits = digits), .groups = "drop") %>%
    mutate(
      section = "Organ-score components",
      characteristic = characteristic,
      row_order = row_order
    ) %>%
    select(section, characteristic, row_order, organ, value)
}

component_binary_row <- function(dat, org, var, characteristic, row_order) {
  dat %>%
    filter(WL_ORG == org) %>%
    group_by(organ) %>%
    summarise(value = fmt_n_pct(sum(.data[[var]] == 1L, na.rm = TRUE), n()), .groups = "drop") %>%
    mutate(
      section = "Organ-score components",
      characteristic = characteristic,
      row_order = row_order
    ) %>%
    select(section, characteristic, row_order, organ, value)
}

recode_sex <- function(x) {
  case_when(
    x == "F" ~ "Female",
    x == "M" ~ "Male",
    is.na(x) | x == "" ~ "Unknown",
    TRUE ~ "Unknown"
  )
}

recode_race <- function(x) {
  race_upper <- str_to_upper(coalesce(x, ""))
  case_when(
    str_detect(race_upper, "WHITE") ~ "White",
    str_detect(race_upper, "BLACK") ~ "Black",
    str_detect(race_upper, "ASIAN") ~ "Asian",
    TRUE ~ "Other/Unknown"
  )
}

message("Reading primary waitlist-period pollution analysis dataset")
analysis_dat <- read_csv(in_path, show_col_types = FALSE) %>%
  filter(WL_ORG %in% names(organ_labels)) %>%
  mutate(
    organ = factor(recode(WL_ORG, !!!organ_labels), levels = organ_levels),
    sex_group = factor(recode_sex(sex), levels = sex_levels),
    race_group = factor(recode_race(race), levels = race_levels),
    listing_year = as.integer(listing_year),
    organ_score_display = case_when(
      WL_ORG == "LI" ~ organ_score - 6200,
      TRUE ~ organ_score
    )
  )

denoms <- analysis_dat %>%
  count(organ, name = "n")

scalar_rows <- bind_rows(
  analysis_dat %>%
    group_by(organ) %>%
    summarise(value = fmt_int(n()), .groups = "drop") %>%
    mutate(section = "Cohort", characteristic = "Waitlist listings, n", row_order = 1),
  analysis_dat %>%
    group_by(organ) %>%
    summarise(value = fmt_int(n_distinct(PERS_ID)), .groups = "drop") %>%
    mutate(section = "Cohort", characteristic = "Unique candidates, n", row_order = 2),
  analysis_dat %>%
    group_by(organ) %>%
    summarise(value = fmt_int(n_distinct(listing_center)), .groups = "drop") %>%
    mutate(section = "Cohort", characteristic = "Listing centers, n", row_order = 3),
  analysis_dat %>%
    group_by(organ) %>%
    summarise(value = fmt_median_iqr(age, digits = 1), .groups = "drop") %>%
    mutate(section = "Baseline characteristics", characteristic = "Age at listing, median (IQR), years", row_order = 4),
  analysis_dat %>%
    group_by(organ) %>%
    summarise(value = fmt_median_iqr(listing_year, digits = 0), .groups = "drop") %>%
    mutate(section = "Baseline characteristics", characteristic = "Listing year, median (IQR)", row_order = 5)
) %>%
  select(section, characteristic, row_order, organ, value)

score_rows <- bind_rows(
  analysis_dat %>%
    filter(WL_ORG == "HR") %>%
    group_by(organ) %>%
    summarise(value = fmt_median_iqr(organ_score_display, digits = 1), .groups = "drop") %>%
    mutate(section = "Organ-specific baseline severity", characteristic = "Heart US-CRS proxy, median (IQR)", row_order = 6),
  analysis_dat %>%
    filter(WL_ORG == "LI") %>%
    group_by(organ) %>%
    summarise(value = fmt_median_iqr(organ_score_display, digits = 0), .groups = "drop") %>%
    mutate(section = "Organ-specific baseline severity", characteristic = "Liver MELD/PELD score, median (IQR)", row_order = 8),
  analysis_dat %>%
    filter(WL_ORG == "LU") %>%
    group_by(organ) %>%
    summarise(value = fmt_median_iqr(organ_score_display, digits = 1), .groups = "drop") %>%
    mutate(section = "Organ-specific baseline severity", characteristic = "Lung LAS/CAS component proxy, median (IQR)", row_order = 9)
) %>%
  select(section, characteristic, row_order, organ, value)

score_component_rows <- bind_rows(
  component_cont_row(analysis_dat, "HR", "hr_albumin_imputed", "Heart albumin score input, median (IQR), g/dL", 10, digits = 1),
  component_cont_row(analysis_dat, "HR", "hr_bilirubin_imputed", "Heart bilirubin score input, median (IQR), mg/dL", 11, digits = 1),
  component_cont_row(analysis_dat, "HR", "hr_egfr", "Heart eGFR score input, median (IQR), mL/min/1.73 m2", 12, digits = 1),
  component_cont_row(analysis_dat, "HR", "hr_sodium_imputed", "Heart sodium score input, median (IQR), mEq/L", 13, digits = 1),
  component_binary_row(analysis_dat, "HR", "hr_durable_lvad", "Heart durable LVAD score input, n (%)", 14),
  component_binary_row(analysis_dat, "HR", "hr_short_mcs", "Heart short-term MCS score input, n (%)", 15),
  component_cont_row(analysis_dat, "HR", "hr_bnp_imputed", "Heart BNP score input, median (IQR), pg/mL", 16, digits = 0),
  component_binary_row(analysis_dat, "KI", "kidney_no_dialysis_time", "Kidney no dialysis time covariate, n (%)", 17),
  component_cont_row(
    analysis_dat,
    "KI",
    "kidney_dialysis_years",
    "Kidney dialysis duration covariate among candidates with positive time, median (IQR), years",
    18,
    digits = 1,
    transform = function(x) if_else(x > 0, x, NA_real_)
  ),
  component_binary_row(analysis_dat, "KI", "kidney_diabetes", "Kidney diabetes covariate, n (%)", 19),
  component_cont_row(analysis_dat, "LI", "liver_meld", "Liver MELD/PELD score input, median (IQR)", 20, digits = 0, transform = function(x) x - 6200),
  component_cont_row(analysis_dat, "LU", "lung_fev1_score_input", "Lung FEV1 score input, median (IQR)", 21, digits = 1),
  component_cont_row(analysis_dat, "LU", "lung_fvc_score_input", "Lung FVC score input, median (IQR)", 22, digits = 1),
  component_cont_row(analysis_dat, "LU", "lung_pco2_score_input", "Lung PCO2 score input, median (IQR), mm Hg", 23, digits = 1),
  component_cont_row(analysis_dat, "LU", "lung_resting_o2_score_input", "Lung resting oxygen score input, median (IQR)", 24, digits = 1),
  component_cont_row(analysis_dat, "LU", "lung_six_min_walk_score_input", "Lung 6-minute walk score input, median (IQR), ft", 25, digits = 0),
  component_cont_row(analysis_dat, "LU", "lung_pulm_art_mean_score_input", "Lung pulmonary artery mean pressure score input, median (IQR), mm Hg", 26, digits = 1),
  component_cont_row(analysis_dat, "LU", "lung_cardiac_output_score_input", "Lung cardiac output score input, median (IQR), L/min", 27, digits = 1),
  component_binary_row(analysis_dat, "LU", "lung_ventilator_score_input", "Lung ventilator score input, n (%)", 28),
  component_binary_row(analysis_dat, "LU", "lung_ecmo_score_input", "Lung ECMO score input, n (%)", 29),
  component_binary_row(analysis_dat, "LU", "lung_corticosteroid_score_input", "Lung corticosteroid dependence score input, n (%)", 30)
)

exposure_rows <- bind_rows(
  analysis_dat %>%
    filter(listing_year <= 2023, is.finite(pm25_waitlist_ug_m3)) %>%
    group_by(organ) %>%
    summarise(value = fmt_median_iqr_range(pm25_waitlist_ug_m3, digits = 1), .groups = "drop") %>%
    mutate(section = "Waitlist-period pollution exposure", characteristic = "PM2.5, median (IQR); range, ug/m3", row_order = 80),
  analysis_dat %>%
    filter(listing_year <= 2023, is.finite(o3_waitlist_ppb)) %>%
    group_by(organ) %>%
    summarise(value = fmt_median_iqr_range(o3_waitlist_ppb, digits = 1), .groups = "drop") %>%
    mutate(section = "Waitlist-period pollution exposure", characteristic = "O3, median (IQR); range, ppb", row_order = 81),
  analysis_dat %>%
    filter(listing_year <= 2025, is.finite(no2_waitlist)) %>%
    group_by(organ) %>%
    summarise(value = fmt_median_iqr_range(no2_waitlist, digits = 1), .groups = "drop") %>%
    mutate(section = "Waitlist-period pollution exposure", characteristic = "NO2, median (IQR); range, ppb", row_order = 82)
) %>%
  select(section, characteristic, row_order, organ, value)

sex_rows <- analysis_dat %>%
  count(organ, sex_group, name = "n") %>%
  complete(organ, sex_group = factor(sex_levels, levels = sex_levels), fill = list(n = 0)) %>%
  left_join(denoms, by = "organ", suffix = c("", "_denom")) %>%
  mutate(
    section = "Sex, n (%)",
    characteristic = as.character(sex_group),
    row_order = 40 + as.integer(sex_group),
    value = mapply(fmt_n_pct, n, n_denom)
  ) %>%
  select(section, characteristic, row_order, organ, value)

race_rows <- analysis_dat %>%
  count(organ, race_group, name = "n") %>%
  complete(organ, race_group = factor(race_levels, levels = race_levels), fill = list(n = 0)) %>%
  left_join(denoms, by = "organ", suffix = c("", "_denom")) %>%
  mutate(
    section = "Race, n (%)",
    characteristic = as.character(race_group),
    row_order = 50 + as.integer(race_group),
    value = mapply(fmt_n_pct, n, n_denom)
  ) %>%
  select(section, characteristic, row_order, organ, value)

outcome_rows <- analysis_dat %>%
  mutate(event_type = factor(event_type, levels = outcome_levels)) %>%
  count(organ, event_type, name = "n") %>%
  complete(organ, event_type = factor(outcome_levels, levels = outcome_levels), fill = list(n = 0)) %>%
  left_join(denoms, by = "organ", suffix = c("", "_denom")) %>%
  mutate(
    section = "Waitlist outcome, n (%)",
    characteristic = recode(as.character(event_type), !!!outcome_labels),
    row_order = 60 + as.integer(event_type),
    value = mapply(fmt_n_pct, n, n_denom)
  ) %>%
  select(section, characteristic, row_order, organ, value)

table_long <- bind_rows(scalar_rows, score_rows, score_component_rows, sex_rows, race_rows, outcome_rows, exposure_rows) %>%
  mutate(organ = factor(organ, levels = organ_levels)) %>%
  arrange(row_order, organ)

table_wide <- table_long %>%
  select(section, characteristic, organ, value) %>%
  pivot_wider(names_from = organ, values_from = value) %>%
  arrange(match(section, unique(table_long$section)), match(characteristic, unique(table_long$characteristic))) %>%
  mutate(across(all_of(organ_levels), ~ tidyr::replace_na(.x, "")))

csv_path <- file.path(out_dir, "table1_baseline_characteristics_by_organ.csv")
html_path <- file.path(out_dir, "table1_baseline_characteristics_by_organ.html")
caption_path <- file.path(out_dir, "table1_baseline_characteristics_caption.txt")
write_csv(table_wide, csv_path)

gt_table <- table_wide %>%
  gt(groupname_col = "section", rowname_col = "characteristic") %>%
  tab_header(
    title = "Baseline Characteristics of Waitlist Candidates by Organ",
    subtitle = "Primary waitlist-period pollution analysis cohort"
  ) %>%
  cols_label(
    Heart = "Heart",
    Kidney = "Kidney",
    Liver = "Liver",
    Lung = "Lung"
  ) %>%
  tab_source_note(
    source_note = "Values are median (IQR) for continuous variables and n (%) for categorical variables unless otherwise indicated."
  ) %>%
  tab_source_note(
    source_note = "Organ-specific severity scores and score components are shown only for the organ to which each score applies. Heart score components use median-imputed laboratory values where needed; lung score components reflect the values used in the proxy score, with missing numeric components set to 0. Liver MELD/PELD is displayed after subtracting the SRTR 6200 offset from the stored MELD/PELD field."
  ) %>%
  tab_source_note(
    source_note = "Pollution exposures are day-weighted waitlist-period means. PM2.5 and O3 are summarized among candidates with exposure follow-up through 2023; NO2 is summarized among candidates with exposure follow-up through 2025."
  ) %>%
  tab_source_note(
    source_note = "Waitlist outcomes are mutually exclusive final observed outcomes in the primary cohort before pollutant-specific exposure censoring."
  ) %>%
  tab_options(
    table.font.names = "Arial",
    table.font.size = px(13),
    heading.title.font.size = px(18),
    heading.subtitle.font.size = px(13),
    row_group.font.weight = "bold",
    column_labels.font.weight = "bold",
    table.width = pct(100),
    data_row.padding = px(4)
  )

gtsave(gt_table, html_path)

writeLines(
  c(
    "Table 1 caption:",
    "Baseline characteristics of waitlist listings included in the primary waitlist-period pollution analysis cohort, stratified by listed organ.",
    "Values are median (IQR) for continuous variables and n (%) for categorical variables unless otherwise indicated.",
    "Organ-specific severity scores and score components are shown only for the organ to which each score applies. Heart score components use median-imputed laboratory values where needed; lung score components reflect the values used in the proxy score, with missing numeric components set to 0. Liver MELD/PELD is displayed after subtracting the SRTR 6200 offset from the stored MELD/PELD field.",
    "Waitlist outcomes are mutually exclusive final observed outcomes in the primary cohort before pollutant-specific exposure censoring.",
    "Pollution exposures are day-weighted waitlist-period means; exposure rows include median (IQR) and range."
  ),
  caption_path
)

message("Wrote Table 1 to ", normalizePath(out_dir, winslash = "/"))
