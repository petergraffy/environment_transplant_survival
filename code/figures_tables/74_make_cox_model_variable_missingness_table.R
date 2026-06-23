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

analysis_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_analysis_dataset.csv.gz"
)
out_dir <- file.path("output", "tables", "primary_cox_model_variables_missingness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")
organ_order <- unname(organ_labels)

fmt_int <- function(x) {
  formatC(as.integer(round(x)), format = "d", big.mark = ",")
}

fmt_pct <- function(x) {
  ifelse(is.na(x), "", sprintf("%.1f%%", x))
}

missing_summary <- function(dat, org, variable, label, role, modeling_form, handling, denominator_note, section) {
  org_dat <- dat %>% filter(WL_ORG == org)
  value <- org_dat[[variable]]
  tibble(
    section = section,
    organ = unname(organ_labels[[org]]),
    variable = label,
    role = role,
    modeling_form_or_weighting = modeling_form,
    baseline_missingness_denominator = denominator_note,
    n = nrow(org_dat),
    missing_n = sum(is.na(value)),
    missing_percent = 100 * missing_n / n,
    model_handling = handling
  )
}

exposure_missing_summary <- function(dat, org, variable, label, end_year, modeling_form) {
  org_dat <- dat %>%
    filter(WL_ORG == org, as.integer(listing_year) <= end_year, .data[[paste0("followup_days_", variable)]] > 0)
  exposure_var <- case_when(
    variable == "pm25" ~ "pm25_waitlist_ug_m3",
    variable == "o3" ~ "o3_waitlist_ppb",
    variable == "no2" ~ "no2_waitlist"
  )
  value <- org_dat[[exposure_var]]
  tibble(
    section = "Exposure",
    organ = unname(organ_labels[[org]]),
    variable = label,
    role = "Primary exposure; modeled independently from other pollutants",
    modeling_form_or_weighting = modeling_form,
    baseline_missingness_denominator = paste0("Person waitlist episodes through ", end_year, " with positive pollutant-specific follow-up"),
    n = nrow(org_dat),
    missing_n = sum(is.na(value)),
    missing_percent = 100 * missing_n / n,
    model_handling = "Rows with missing pollutant exposure were excluded from that pollutant-specific Cox model"
  )
}

survival_variable_summary <- function(dat, org, variable, label, pollutant_label, end_year, modeling_form) {
  org_dat <- dat %>%
    filter(WL_ORG == org, as.integer(listing_year) <= end_year, .data[[paste0("followup_days_", pollutant_label)]] > 0)
  value <- org_dat[[variable]]
  tibble(
    section = "Outcome and time scale",
    organ = unname(organ_labels[[org]]),
    variable = label,
    role = "Survival outcome/time variable",
    modeling_form_or_weighting = modeling_form,
    baseline_missingness_denominator = paste0("Person waitlist episodes through ", end_year, " with positive pollutant-specific follow-up"),
    n = nrow(org_dat),
    missing_n = sum(is.na(value)),
    missing_percent = 100 * missing_n / n,
    model_handling = "Required for model fitting"
  )
}

dat <- read_csv(analysis_path, show_col_types = FALSE) %>%
  filter(WL_ORG %in% names(organ_labels)) %>%
  mutate(
    listing_year = as.integer(listing_year),
    sex = na_if(as.character(sex), "Missing"),
    race = na_if(as.character(race), "Missing"),
    listing_center = na_if(as.character(listing_center), "Missing")
  )

common_rows <- bind_rows(lapply(names(organ_labels), function(org) {
  bind_rows(
    missing_summary(dat, org, "age", "Age at listing", "Adjustment covariate", "Continuous, linear", "Complete cases required for model fitting", "All primary cohort person waitlist episodes for organ", "Core covariates"),
    missing_summary(dat, org, "sex", "Sex", "Adjustment covariate", "Categorical indicator variables", "Missing values modeled as an explicit Missing category", "All primary cohort person waitlist episodes for organ", "Core covariates"),
    missing_summary(dat, org, "race", "Race", "Adjustment covariate", "Categorical indicator variables", "Missing values modeled as an explicit Missing category", "All primary cohort person waitlist episodes for organ", "Core covariates"),
    missing_summary(dat, org, "listing_year", "Listing year", "Temporal adjustment covariate", "Categorical indicator variables", "Complete cases required for model fitting", "All primary cohort person waitlist episodes for organ", "Core covariates"),
    missing_summary(dat, org, "listing_center", "Listing transplant center", "Center adjustment", "Stratified baseline hazard using strata(listing_center)", "Missing values modeled as an explicit Missing stratum", "All primary cohort person waitlist episodes for organ", "Core covariates")
  )
}))

survival_rows <- bind_rows(lapply(names(organ_labels), function(org) {
  bind_rows(
    survival_variable_summary(dat, org, "followup_days_pm25", "PM2.5 model follow-up time", "pm25", 2023L, "Days from waitlist start to adverse event, competing/nonadverse waitlist exit, administrative censoring, or PM2.5 data end"),
    survival_variable_summary(dat, org, "event_pm25", "PM2.5 model adverse-event indicator", "pm25", 2023L, "Death or delisting due to deterioration before PM2.5 data end"),
    survival_variable_summary(dat, org, "followup_days_o3", "O3 model follow-up time", "o3", 2023L, "Days from waitlist start to adverse event, competing/nonadverse waitlist exit, administrative censoring, or O3 data end"),
    survival_variable_summary(dat, org, "event_o3", "O3 model adverse-event indicator", "o3", 2023L, "Death or delisting due to deterioration before O3 data end"),
    survival_variable_summary(dat, org, "followup_days_no2", "NO2 model follow-up time", "no2", 2025L, "Days from waitlist start to adverse event, competing/nonadverse waitlist exit, administrative censoring, or NO2 data end"),
    survival_variable_summary(dat, org, "event_no2", "NO2 model adverse-event indicator", "no2", 2025L, "Death or delisting due to deterioration before NO2 data end")
  )
}))

exposure_rows <- bind_rows(lapply(names(organ_labels), function(org) {
  bind_rows(
    exposure_missing_summary(dat, org, "pm25", "PM2.5 waitlist-period exposure", 2023L, "Day-weighted mean exposure during observed waitlist follow-up; HR per 5-ug/m3 increase"),
    exposure_missing_summary(dat, org, "o3", "O3 waitlist-period exposure", 2023L, "Day-weighted mean exposure during observed waitlist follow-up; HR per 10-ppb increase"),
    exposure_missing_summary(dat, org, "no2", "NO2 waitlist-period exposure", 2025L, "Day-weighted mean exposure during observed waitlist follow-up; HR per 10-ppb increase")
  )
}))

organ_score_rows <- bind_rows(
  missing_summary(dat, "HR", "organ_score", "Heart US-CRS proxy", "Baseline severity adjustment", "Continuous, linear composite score", "Derived from available baseline fields after median imputation of laboratory components", "Heart person waitlist episodes", "Organ-specific baseline severity"),
  missing_summary(dat, "KI", "kidney_no_dialysis_time", "No dialysis time", "Kidney baseline severity adjustment", "Binary indicator", "Constructed as 1 for no dialysis, missing dialysis date, or 0 dialysis duration", "Kidney person waitlist episodes", "Organ-specific baseline severity"),
  missing_summary(dat, "KI", "kidney_dialysis_years", "Dialysis duration", "Kidney baseline severity adjustment", "Continuous, linear years", "Set to 0 when no dialysis time indicator was 1", "Kidney person waitlist episodes", "Organ-specific baseline severity"),
  missing_summary(dat, "KI", "kidney_diabetes", "Diabetes", "Kidney baseline severity adjustment", "Binary indicator", "Constructed from SRTR diabetes fields", "Kidney person waitlist episodes", "Organ-specific baseline severity"),
  missing_summary(dat, "LI", "organ_score", "Liver MELD/PELD score", "Baseline severity adjustment", "Continuous, linear score", "Initial MELD/PELD used, with last MELD/PELD used when initial value was unavailable", "Liver person waitlist episodes", "Organ-specific baseline severity"),
  missing_summary(dat, "LU", "organ_score", "Lung LAS/CAS component proxy", "Baseline severity adjustment", "Continuous, linear composite score", "Derived from baseline lung components; missing numeric components set to 0 in proxy score", "Lung person waitlist episodes", "Organ-specific baseline severity")
)

component_rows <- bind_rows(
  missing_summary(dat, "HR", "hr_albumin_observed", "Heart albumin", "Heart US-CRS proxy component", "Median-imputed continuous input to proxy score", "Raw missingness shown; missing values median-imputed within heart cohort for score calculation", "Heart person waitlist episodes", "Organ-score components"),
  missing_summary(dat, "HR", "hr_bilirubin_observed", "Heart bilirubin", "Heart US-CRS proxy component", "Median-imputed log1p continuous input to proxy score", "Raw missingness shown; missing values median-imputed within heart cohort for score calculation", "Heart person waitlist episodes", "Organ-score components"),
  missing_summary(dat, "HR", "hr_sodium_observed", "Heart sodium", "Heart US-CRS proxy component", "Median-imputed continuous input to proxy score", "Raw missingness shown; missing values median-imputed within heart cohort for score calculation", "Heart person waitlist episodes", "Organ-score components"),
  missing_summary(dat, "HR", "hr_bnp_observed", "Heart BNP", "Heart US-CRS proxy component", "Median-imputed log1p continuous input to proxy score", "Raw missingness shown; missing values median-imputed within heart cohort for score calculation", "Heart person waitlist episodes", "Organ-score components"),
  missing_summary(dat, "HR", "hr_durable_lvad", "Durable LVAD", "Heart US-CRS proxy component", "Binary indicator", "Constructed from CAND_THOR and candidate VAD brand codes; missing/no qualifying device coded 0", "Heart person waitlist episodes", "Organ-score components"),
  missing_summary(dat, "HR", "hr_short_mcs", "Short-term mechanical circulatory support", "Heart US-CRS proxy component", "Binary indicator", "Constructed from ECMO/RVAD fields; missing/no support coded 0", "Heart person waitlist episodes", "Organ-score components"),
  missing_summary(dat, "HR", "hr_egfr", "Heart eGFR", "Heart US-CRS proxy component", "Continuous input derived with 2021 CKD-EPI creatinine equation", "Creatinine was median-imputed before eGFR calculation when missing", "Heart person waitlist episodes", "Organ-score components"),
  missing_summary(dat, "LU", "lung_fev1_observed", "Lung FEV1", "Lung LAS/CAS proxy component", "Continuous input to proxy score", "Raw missingness shown; missing values set to 0 for proxy score calculation", "Lung person waitlist episodes", "Organ-score components"),
  missing_summary(dat, "LU", "lung_fvc_observed", "Lung FVC", "Lung LAS/CAS proxy component", "Continuous input to proxy score", "Raw missingness shown; missing values set to 0 for proxy score calculation", "Lung person waitlist episodes", "Organ-score components"),
  missing_summary(dat, "LU", "lung_pco2_observed", "Lung PCO2", "Lung LAS/CAS proxy component", "Continuous input to proxy score", "Raw missingness shown; missing values set to 0 for proxy score calculation", "Lung person waitlist episodes", "Organ-score components"),
  missing_summary(dat, "LU", "lung_resting_o2_observed", "Lung resting oxygen", "Lung LAS/CAS proxy component", "Continuous input to proxy score", "Raw missingness shown; missing values set to 0 for proxy score calculation", "Lung person waitlist episodes", "Organ-score components"),
  missing_summary(dat, "LU", "lung_six_min_walk_observed", "Lung 6-minute walk distance", "Lung LAS/CAS proxy component", "Continuous input to proxy score", "Raw missingness shown; missing values set to 0 for proxy score calculation", "Lung person waitlist episodes", "Organ-score components"),
  missing_summary(dat, "LU", "lung_pulm_art_mean_observed", "Lung mean pulmonary artery pressure", "Lung LAS/CAS proxy component", "Continuous input to proxy score", "Raw missingness shown; missing values set to 0 for proxy score calculation", "Lung person waitlist episodes", "Organ-score components"),
  missing_summary(dat, "LU", "lung_cardiac_output_observed", "Lung cardiac output", "Lung LAS/CAS proxy component", "Continuous input to proxy score", "Raw missingness shown; missing values set to 0 for proxy score calculation", "Lung person waitlist episodes", "Organ-score components"),
  missing_summary(dat, "LU", "lung_ventilator_score_input", "Ventilator use", "Lung LAS/CAS proxy component", "Binary indicator", "Missing/no ventilator coded 0", "Lung person waitlist episodes", "Organ-score components"),
  missing_summary(dat, "LU", "lung_ecmo_score_input", "ECMO use", "Lung LAS/CAS proxy component", "Binary indicator", "Missing/no ECMO coded 0", "Lung person waitlist episodes", "Organ-score components"),
  missing_summary(dat, "LU", "lung_corticosteroid_score_input", "Corticosteroid dependence", "Lung LAS/CAS proxy component", "Binary indicator", "Missing/no corticosteroid dependence coded 0", "Lung person waitlist episodes", "Organ-score components")
)

section_levels <- c(
  "Exposure",
  "Outcome and time scale",
  "Core covariates",
  "Organ-specific baseline severity",
  "Organ-score components"
)

table_long <- bind_rows(exposure_rows, survival_rows, common_rows, organ_score_rows, component_rows) %>%
  mutate(
    section = factor(section, levels = section_levels),
    organ = factor(organ, levels = organ_order),
    missing = paste0(fmt_int(missing_n), "/", fmt_int(n), " (", fmt_pct(missing_percent), ")")
  ) %>%
  arrange(section, organ, variable)

write_csv(table_long, file.path(out_dir, "cox_model_variables_missingness_long.csv"))

table_display <- table_long %>%
  transmute(
    Section = section,
    Organ = organ,
    Variable = variable,
    Role = role,
    `Modeling form / weighting` = modeling_form_or_weighting,
    `Missing, n/N (%)` = missing,
    `Missingness denominator` = baseline_missingness_denominator,
    `Model handling` = model_handling
  )

write_csv(table_display, file.path(out_dir, "cox_model_variables_missingness_supplement_table.csv"))

gt_table <- table_display %>%
  gt(groupname_col = "Section") %>%
  tab_header(title = "Variables Included in Primary Cause-Specific Cox Models") %>%
  tab_source_note(source_note = "Missingness is reported in the primary baseline analysis cohort for each organ unless otherwise noted. Exposure missingness is pollutant-specific and uses the cohort eligible for that pollutant's data period.") %>%
  tab_source_note(source_note = "PM2.5 and O3 exposure data were available through 2023; NO2 exposure data were available through 2025.") %>%
  tab_options(
    table.font.names = "Arial",
    table.font.size = px(11),
    heading.title.font.size = px(16),
    row_group.font.weight = "bold",
    column_labels.font.weight = "bold",
    data_row.padding = px(3),
    table.width = pct(100)
  )

gtsave(gt_table, file.path(out_dir, "cox_model_variables_missingness_supplement_table.html"))

writeLines(
  c(
    "Supplemental table note:",
    "This table describes variables included in the primary cause-specific Cox models. Each pollutant was modeled in a separate organ-specific Cox model. PM2.5 and O3 were modeled using day-weighted waitlist-period exposure through 2023; NO2 was modeled using day-weighted waitlist-period exposure through 2025. Missingness is shown before model-specific complete-case exclusion or imputation/derived-variable handling, except for constructed binary indicators where missing/no qualifying feature was coded as 0."
  ),
  file.path(out_dir, "cox_model_variables_missingness_table_note.txt")
)

message("Wrote Cox model variable missingness table to ", normalizePath(out_dir, winslash = "/"))
