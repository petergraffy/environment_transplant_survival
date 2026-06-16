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
    filter(WL_ORG == "KI") %>%
    group_by(organ) %>%
    summarise(value = fmt_median_iqr(organ_score_display, digits = 1), .groups = "drop") %>%
    mutate(section = "Organ-specific baseline severity", characteristic = "Kidney dialysis-age score, median (IQR)", row_order = 7),
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

exposure_rows <- bind_rows(
  analysis_dat %>%
    filter(listing_year <= 2023, is.finite(pm25_waitlist_ug_m3)) %>%
    group_by(organ) %>%
    summarise(value = fmt_median_iqr_range(pm25_waitlist_ug_m3, digits = 1), .groups = "drop") %>%
    mutate(section = "Waitlist-period pollution exposure", characteristic = "PM2.5, median (IQR); range, ug/m3", row_order = 30),
  analysis_dat %>%
    filter(listing_year <= 2023, is.finite(o3_waitlist_ppb)) %>%
    group_by(organ) %>%
    summarise(value = fmt_median_iqr_range(o3_waitlist_ppb, digits = 1), .groups = "drop") %>%
    mutate(section = "Waitlist-period pollution exposure", characteristic = "O3, median (IQR); range, ppb", row_order = 31),
  analysis_dat %>%
    filter(listing_year <= 2025, is.finite(no2_waitlist)) %>%
    group_by(organ) %>%
    summarise(value = fmt_median_iqr_range(no2_waitlist, digits = 1), .groups = "drop") %>%
    mutate(section = "Waitlist-period pollution exposure", characteristic = "NO2, median (IQR); range, ppb", row_order = 32)
) %>%
  select(section, characteristic, row_order, organ, value)

sex_rows <- analysis_dat %>%
  count(organ, sex_group, name = "n") %>%
  complete(organ, sex_group = factor(sex_levels, levels = sex_levels), fill = list(n = 0)) %>%
  left_join(denoms, by = "organ", suffix = c("", "_denom")) %>%
  mutate(
    section = "Sex, n (%)",
    characteristic = as.character(sex_group),
    row_order = 10 + as.integer(sex_group),
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
    row_order = 20 + as.integer(race_group),
    value = mapply(fmt_n_pct, n, n_denom)
  ) %>%
  select(section, characteristic, row_order, organ, value)

table_long <- bind_rows(scalar_rows, score_rows, sex_rows, race_rows, exposure_rows) %>%
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
    source_note = "Organ-specific severity scores are shown only for the organ to which each score applies. Liver MELD/PELD is displayed after subtracting the SRTR 6200 offset from the stored MELD/PELD field."
  ) %>%
  tab_source_note(
    source_note = "Pollution exposures are day-weighted waitlist-period means. PM2.5 and O3 are summarized among candidates with exposure follow-up through 2023; NO2 is summarized among candidates with exposure follow-up through 2025."
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
    "Organ-specific severity scores are shown only for the organ to which each score applies. Liver MELD/PELD is displayed after subtracting the SRTR 6200 offset from the stored MELD/PELD field.",
    "Pollution exposures are day-weighted waitlist-period means; exposure rows include median (IQR) and range."
  ),
  caption_path
)

message("Wrote Table 1 to ", normalizePath(out_dir, winslash = "/"))
