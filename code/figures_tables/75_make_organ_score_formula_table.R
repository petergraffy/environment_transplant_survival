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
  library(tibble)
})

out_dir <- file.path("output", "tables", "organ_score_formula_weights")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

formula_table <- tribble(
  ~Organ, ~`Model variable`, ~`Component`, ~`Source variable(s)`, ~`Transformation`, ~`Coefficient / weight`, ~`Missing-data handling used for score`,
  "Heart", "Heart US-CRS proxy", "Albumin", "RiskStratDataHR HrSevFailAlbumin; fallback statjust_hr1a CANHX_LAB_ALBUMIN; fallback cand_thor CAN_TOT_ALBUMIN", "Continuous, g/dL", "-0.656", "Median imputed within heart cohort before score calculation",
  "Heart", "Heart US-CRS proxy", "Bilirubin", "RiskStratDataHR HrSevFailBilirubin; fallback statjust_hr1a CANHX_LAB_BILI; fallback CAN_TOT_BILI when available", "log1p(value)", "0.617", "Median imputed within heart cohort before score calculation",
  "Heart", "Heart US-CRS proxy", "eGFR", "RiskStratDataHR HrSevFailCreatinine; fallback statjust_hr1a CANHX_LAB_SERUM_CREAT; fallback cand_thor CAN_MOST_RECENT_CREAT", "2021 CKD-EPI creatinine eGFR", "-0.012", "Creatinine median imputed within heart cohort before eGFR calculation",
  "Heart", "Heart US-CRS proxy", "Sodium", "RiskStratDataHR HrSevFailSodium; fallback statjust_hr1a CANHX_LAB_SODIUM; fallback CAN_LAST_SERUM_SODIUM when available", "Continuous, mEq/L", "-0.077", "Median imputed within heart cohort before score calculation",
  "Heart", "Heart US-CRS proxy", "Durable LVAD", "CAND_THOR CAN_VAD1 and CAN_VAD2 durable device brand codes", "Binary indicator", "-0.377", "Missing/no qualifying durable device coded 0",
  "Heart", "Heart US-CRS proxy", "Short-term mechanical circulatory support", "statjust_hr1a CANHX_ECMO or CANHX_RVAD_TYPE; fallback cand_thor CAN_ECMO", "Binary indicator", "1.092", "Missing/no support coded 0",
  "Heart", "Heart US-CRS proxy", "BNP", "RiskStratDataHR HrSevFailBnp; fallback statjust_hr1a CANHX_LAB_BNP", "log1p(value)", "0.433", "Median imputed within heart cohort before score calculation",
  "Kidney", "No composite organ score", "No dialysis time", "CAN_ON_DIAL/CAN_DIAL and CAN_DIAL_DT", "Binary indicator", "Estimated model coefficient", "Constructed as 1 for not on dialysis, missing dialysis date, or 0 dialysis duration",
  "Kidney", "No composite organ score", "Dialysis duration", "CAN_DIAL_DT and waitlist episode start date", "Continuous years", "Estimated model coefficient per year", "Set to 0 when no dialysis time indicator was 1",
  "Kidney", "No composite organ score", "Diabetes", "CAN_DIAB or CAN_DIAB_TY", "Binary indicator", "Estimated model coefficient", "Constructed from SRTR diabetes fields",
  "Liver", "MELD/PELD score", "Initial MELD/PELD", "CAN_INIT_SRTR_LAB_MELD; fallback CAN_LAST_SRTR_LAB_MELD", "Continuous score as stored in SRTR model; Table 1 display subtracts 6200 offset", "Estimated model coefficient per 1-point increase", "Initial value used when available; last value used when initial missing",
  "Lung", "Lung LAS/CAS component proxy", "FEV1", "CAN_FEV1", "Continuous", "-0.010", "Missing numeric values set to 0 for proxy score calculation",
  "Lung", "Lung LAS/CAS component proxy", "FVC", "CAN_FVC", "Continuous", "-0.005", "Missing numeric values set to 0 for proxy score calculation",
  "Lung", "Lung LAS/CAS component proxy", "PCO2", "CAN_PCO2", "Continuous, mm Hg", "0.020", "Missing numeric values set to 0 for proxy score calculation",
  "Lung", "Lung LAS/CAS component proxy", "Resting oxygen requirement", "CAN_AT_REST_O2", "Continuous", "0.010", "Missing numeric values set to 0 for proxy score calculation",
  "Lung", "Lung LAS/CAS component proxy", "6-minute walk distance", "CAN_SIX_MIN_WALK", "Continuous, ft", "-0.002", "Missing numeric values set to 0 for proxy score calculation",
  "Lung", "Lung LAS/CAS component proxy", "Mean pulmonary artery pressure", "CAN_PULM_ART_MEAN", "Continuous, mm Hg", "0.010", "Missing numeric values set to 0 for proxy score calculation",
  "Lung", "Lung LAS/CAS component proxy", "Cardiac output", "CAN_CARDIAC_OUTPUT", "Continuous, L/min", "-0.050", "Missing numeric values set to 0 for proxy score calculation",
  "Lung", "Lung LAS/CAS component proxy", "Ventilator use", "CAN_VENTILATOR or CAN_ON_VENTILATOR", "Binary indicator", "1.000", "Missing/no ventilator coded 0",
  "Lung", "Lung LAS/CAS component proxy", "ECMO use", "CAN_ECMO", "Binary indicator", "1.500", "Missing/no ECMO coded 0",
  "Lung", "Lung LAS/CAS component proxy", "Corticosteroid dependence", "CAN_CORTICOST_DEPND", "Binary indicator", "0.300", "Missing/no corticosteroid dependence coded 0"
)

formula_notes <- tribble(
  ~Organ, ~`Formula / model use`,
  "Heart", "Heart US-CRS proxy = -0.656*albumin + 0.617*log1p(bilirubin) - 0.012*eGFR - 0.077*sodium - 0.377*durable LVAD + 1.092*short-term mechanical circulatory support + 0.433*log1p(BNP).",
  "Kidney", "Kidney models did not use a composite organ score. Instead, no dialysis time, dialysis duration in years, and diabetes were included as separate adjustment covariates.",
  "Liver", "Liver models used the SRTR MELD/PELD field as a continuous adjustment covariate, using initial MELD/PELD with last MELD/PELD as fallback.",
  "Lung", "Lung LAS/CAS component proxy = -0.010*FEV1 - 0.005*FVC + 0.020*PCO2 + 0.010*resting oxygen - 0.002*6-minute walk + 0.010*mean pulmonary artery pressure - 0.050*cardiac output + 1.000*ventilator + 1.500*ECMO + 0.300*corticosteroid dependence."
)

write_csv(formula_table, file.path(out_dir, "organ_score_formula_weights.csv"))
write_csv(formula_notes, file.path(out_dir, "organ_score_formula_notes.csv"))

formula_gt <- formula_table %>%
  gt(groupname_col = "Organ") %>%
  tab_header(title = "Organ-Specific Baseline Severity Scores and Component Weights") %>%
  tab_source_note(source_note = "Heart and lung scores are proxy severity scores derived from available SRTR fields. Kidney models used separate clinical covariates rather than a composite score.") %>%
  tab_source_note(source_note = "log1p(value) denotes log(value + 1).") %>%
  tab_options(
    table.font.names = "Arial",
    table.font.size = px(11),
    heading.title.font.size = px(16),
    row_group.font.weight = "bold",
    column_labels.font.weight = "bold",
    data_row.padding = px(3),
    table.width = pct(100)
  )

notes_gt <- formula_notes %>%
  gt() %>%
  tab_header(title = "Organ-Specific Baseline Severity Formula Summary") %>%
  tab_options(
    table.font.names = "Arial",
    table.font.size = px(12),
    heading.title.font.size = px(16),
    column_labels.font.weight = "bold",
    data_row.padding = px(4),
    table.width = pct(100)
  )

gtsave(formula_gt, file.path(out_dir, "organ_score_formula_weights.html"))
gtsave(notes_gt, file.path(out_dir, "organ_score_formula_notes.html"))

message("Wrote organ score formula tables to ", normalizePath(out_dir, winslash = "/"))
