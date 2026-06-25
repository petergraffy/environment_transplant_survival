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
  library(haven)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
})

saf_paths_source <- file.path("code", "saf_paths.R")
if (!file.exists(saf_paths_source)) saf_paths_source <- "saf_paths.R"
source(saf_paths_source)
saf_paths <- get_saf_paths(release = "q1_2026")

pubsaf_dir <- saf_paths$pubsaf_dir
out_dir <- file.path("output", "figures", "cohort_flow_diagram")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

collapse_path <- file.path(
  "output", "primary_waitlist_period_pollution_cox",
  "overlapping_multilist_collapse_summary.csv"
)
results_path <- file.path(
  "output", "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_cox_results.csv"
)
analysis_dataset_path <- file.path(
  "output", "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_analysis_dataset.csv.gz"
)

collapse_summary <- read_csv(collapse_path, show_col_types = FALSE)
cox_results <- read_csv(results_path, show_col_types = FALSE)

fmt_n <- function(x) comma(as.numeric(x), accuracy = 1)

clean_zip <- function(x) {
  z <- str_extract(as.character(x), "[0-9]+")
  z <- if_else(nchar(z) > 5L, substr(z, 1, 5), z)
  if_else(is.na(z), NA_character_, str_pad(z, 5, side = "left", pad = "0"))
}

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

cohort_accounting_path <- file.path(out_dir, "cohort_flow_accounting.csv")

log_msg("Computing pre-ZCTA cohort accounting from SAF")
candidate_files <- tribble(
  ~candidate_group, ~path,
  "kidney_pancreas", file.path(pubsaf_dir, "cand_kipa.sas7bdat"),
  "liver_intestine", file.path(pubsaf_dir, "cand_liin.sas7bdat"),
  "thoracic", file.path(pubsaf_dir, "cand_thor.sas7bdat")
)

candidate_cols <- c(
  "PERS_ID", "PX_ID", "WL_ORG", "CAN_LISTING_DT", "CAN_SOURCE",
  "CAN_REM_DT", "CAN_DEATH_DT", "CAN_ENDWLFU"
)

candidate_zip <- read_sas(saf_paths$canzip_file) %>%
  transmute(PERS_ID, PX_ID, candidate_zip = clean_zip(CAN_PERM_ZIP)) %>%
  distinct(PERS_ID, PX_ID, .keep_all = TRUE)

candidate_data <- candidate_files %>%
  rowwise() %>%
  mutate(data = list(read_sas(path, col_select = any_of(candidate_cols)))) %>%
  ungroup() %>%
  select(data) %>%
  tidyr::unnest(data)

max_end_date <- as.Date("2025-12-31")
pre_zcta <- candidate_data %>%
  filter(WL_ORG %in% c("HR", "KI", "LI", "LU")) %>%
  left_join(candidate_zip, by = c("PERS_ID", "PX_ID")) %>%
  mutate(
    index_date = CAN_LISTING_DT,
    index_year = as.integer(format(index_date, "%Y")),
    raw_event_date = coalesce(CAN_REM_DT, CAN_DEATH_DT, CAN_ENDWLFU),
    observed_end_date = if_else(is.na(raw_event_date) | raw_event_date > max_end_date, max_end_date, raw_event_date)
  ) %>%
  filter(
    index_year >= 2005L,
    index_year <= 2025L,
    CAN_SOURCE %in% c("A", "R"),
    !is.na(index_date),
    !is.na(observed_end_date),
    observed_end_date >= index_date
  )

cohort_accounting <- bind_rows(
  pre_zcta %>%
    summarise(step = "Eligible registrations before ZCTA exclusion", n = n(), .groups = "drop"),
  pre_zcta %>%
    summarise(step = "Excluded: missing residential ZCTA", n = sum(is.na(candidate_zip)), .groups = "drop"),
  pre_zcta %>%
    filter(!is.na(candidate_zip)) %>%
    summarise(step = "Included: residential ZCTA available", n = n(), .groups = "drop")
)
write_csv(cohort_accounting, cohort_accounting_path)

pre_zcta_n <- cohort_accounting %>%
  filter(step == "Eligible registrations before ZCTA exclusion") %>%
  pull(n)
missing_zcta_n <- cohort_accounting %>%
  filter(step == "Excluded: missing residential ZCTA") %>%
  pull(n)
zcta_available_n <- cohort_accounting %>%
  filter(step == "Included: residential ZCTA available") %>%
  pull(n)

total_registrations <- sum(collapse_summary$registration_rows_before_collapse)
total_episodes <- sum(collapse_summary$person_waitlist_episodes_after_collapse)
total_removed <- sum(collapse_summary$removed_duplicate_registration_rows)
multicenter_episodes <- sum(collapse_summary$multicenter_collapsed_episode_count)
multicenter_rows <- sum(collapse_summary$multicenter_collapsed_registration_rows)

organ_counts <- collapse_summary %>%
  transmute(
    organ_label,
    label = paste0(
      organ_label, "\n",
      "Candidates: ", fmt_n(person_waitlist_episodes_after_collapse)
    )
  )

organ_outcomes_path <- file.path(out_dir, "cohort_flow_organ_outcomes.csv")

log_msg("Computing organ-level outcomes from deduplicated analysis dataset")
organ_outcomes <- read_csv(
  analysis_dataset_path,
  col_select = c(WL_ORG, index_date, observed_end_date, adverse_event, transplant_or_improvement),
  show_col_types = FALSE
) %>%
  mutate(
    organ_label = recode(WL_ORG, !!!setNames(c("Heart", "Kidney", "Liver", "Lung"), c("HR", "KI", "LI", "LU"))),
    index_date = as.Date(index_date),
    observed_end_date = as.Date(observed_end_date),
    waitlist_days = as.numeric(observed_end_date - index_date)
  ) %>%
  group_by(organ_label) %>%
  summarise(
    candidates = n(),
    adverse_events = sum(adverse_event == 1L, na.rm = TRUE),
    transplant_or_improvement = sum(transplant_or_improvement == 1L, na.rm = TRUE),
    median_waitlist_days = median(waitlist_days, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(organ_outcomes, organ_outcomes_path)

organ_outcome_counts <- organ_outcomes %>%
  mutate(
    label = paste0(
      "Adverse events:\n",
      fmt_n(adverse_events), "\n",
      "Transplant/improvement:\n",
      fmt_n(transplant_or_improvement), "\n",
      "Median waitlist:\n",
      fmt_n(median_waitlist_days), " days"
    )
  )

duplicate_breakdown <- collapse_summary %>%
  transmute(label = paste0(organ_label, ": ", fmt_n(removed_duplicate_registration_rows))) %>%
  pull(label) %>%
  str_c(collapse = "\n")

box_df <- tibble(
  id = c(
    "eligible", "zcta", "episodes",
    paste0("organ_", organ_counts$organ_label),
    paste0("outcome_", organ_outcome_counts$organ_label)
  ),
  x = c(5.4, 5.4, 5.4, 1.25, 4.0, 6.8, 9.55, 1.25, 4.0, 6.8, 9.55),
  y = c(11.0, 8.75, 6.35, 4.35, 4.35, 4.35, 4.35, 1.45, 1.45, 1.45, 1.45),
  width = c(7.7, 7.7, 7.7, rep(2.35, 4), rep(2.35, 4)),
  height = c(1.25, 1.15, 1.15, rep(1.1, 4), rep(1.85, 4)),
  label = c(
    paste0(
      "SRTR waitlist registrations assessed for eligibility\n",
      "Heart, kidney, liver, or lung; 2005-2025\n",
      "N = ", fmt_n(pre_zcta_n)
    ),
    paste0(
      "Included registrations with residential ZCTA\n",
      "N = ", fmt_n(zcta_available_n)
    ),
    paste0(
      "Deduplicated candidate-organ episodes\n",
      "N = ", fmt_n(total_episodes)
    ),
    organ_counts$label,
    organ_outcome_counts$label
  )
) %>%
  mutate(
    xmin = x - width / 2,
    xmax = x + width / 2,
    ymin = y - height / 2,
    ymax = y + height / 2
  )

exclusion_df <- tibble(
  id = c("missing_zcta", "duplicate_registrations"),
  x = c(10.95, 10.95),
  y = c(8.75, 7.2),
  width = c(3.0, 3.0),
  height = c(1.15, 2.25),
  label = c(
    paste0(
      "Excluded\n",
      "Missing residential ZCTA\n",
      "n = ", fmt_n(missing_zcta_n)
    ),
    paste0(
      "Excluded\n",
      "Duplicate same-candidate,\n",
      "same-organ registration rows\n",
      "n = ", fmt_n(total_removed), "\n",
      duplicate_breakdown
    )
  )
) %>%
  mutate(
    xmin = x - width / 2,
    xmax = x + width / 2,
    ymin = y - height / 2,
    ymax = y + height / 2
  )

segment_df <- tibble(
  x = c(5.4, 5.4, 5.4, 5.4, 5.4, 5.4),
  y = c(10.38, 8.18, 5.78, 5.78, 5.78, 5.78),
  xend = c(5.4, 5.4, 1.25, 4.0, 6.8, 9.55),
  yend = c(9.33, 6.93, 4.9, 4.9, 4.9, 4.9)
)

exclusion_segment_df <- tibble(
  x = c(5.4, 5.4),
  y = c(8.75, 7.2),
  xend = c(9.45, 9.45),
  yend = c(8.75, 7.2)
)

split_segment_df <- tibble(
  x = c(1.25, 4.0, 6.8, 9.55),
  y = rep(3.8, 4),
  xend = c(1.25, 4.0, 6.8, 9.55),
  yend = rep(2.38, 4)
)

flow_plot <- ggplot() +
  geom_segment(
    data = exclusion_segment_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5,
    lineend = "round",
    color = "black"
  ) +
  geom_segment(
    data = segment_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5,
    lineend = "round",
    color = "black"
  ) +
  geom_segment(
    data = split_segment_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5,
    lineend = "round",
    color = "black"
  ) +
  geom_rect(
    data = box_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white",
    color = "black",
    linewidth = 0.65
  ) +
  geom_rect(
    data = exclusion_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#F2F2F2",
    color = "black",
    linewidth = 0.65
  ) +
  geom_text(
    data = box_df,
    aes(x = x, y = y, label = label),
    size = 3.95,
    lineheight = 0.92,
    color = "black"
  ) +
  geom_text(
    data = exclusion_df,
    aes(x = x, y = y, label = label),
    size = 3.7,
    lineheight = 0.92,
    color = "black"
  ) +
  coord_cartesian(xlim = c(0, 12.8), ylim = c(0.15, 11.8), clip = "off") +
  theme_void(base_size = 15) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(18, 18, 18, 18)
  )

ggsave(
  file.path(out_dir, "cohort_flow_diagram.png"),
  flow_plot,
  width = 13,
  height = 11,
  dpi = 600
)

ggsave(
  file.path(out_dir, "cohort_flow_diagram.pdf"),
  flow_plot,
  width = 13,
  height = 11,
  device = cairo_pdf
)

write_csv(
  bind_rows(box_df, exclusion_df) %>% select(id, label),
  file.path(out_dir, "cohort_flow_diagram_labels.csv")
)

message("Wrote cohort flow diagram to ", normalizePath(out_dir, winslash = "/"))
