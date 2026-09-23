source("code/r_runtime.R")
ensure_user_library()
suppressPackageStartupMessages({library(dplyr); library(readr); library(tibble)})
run <- "output/revision_runs/rolling_20260922"
specs <- list(
  baseline = list("prior_year_pollution_cox_svi/prior_year_pollution_cox_svi_results.csv", c("organ", "exposure")),
  baseline_tiers = list("prior_year_initial_listing_cox_model_set/prior_year_initial_listing_cox_model_set_results.csv", c("model_tier", "model_type", "organ", "exposure")),
  baseline_coefficients = list("prior_year_initial_listing_cox_model_set/prior_year_initial_listing_cox_model_set_full_coefficients.csv", c("model_tier", "model_type", "organ", "term")),
  timevarying = list("timevarying_pollution_cox_svi/timevarying_pollution_cox_svi_results.csv", c("organ", "pollutant")),
  subgroups = list("prior_year_and_timevarying_pollution_subgroup_cox/prior_year_and_timevarying_pollution_subgroup_cox_results.csv", c("model", "organ", "pollutant", "subgroup", "subgroup_level")),
  sensitivities = list("updated_methods_sensitivity_models/updated_methods_sensitivity_results.csv", c("analysis", "model", "organ", "exposure")),
  baseline_severity = list("prior_year_pollution_baseline_severity_associations/prior_year_pollution_baseline_severity_associations.csv", c("organ", "outcome", "pollutant")),
  timevarying_severity = list("timevarying_pollution_severity_models_with_listing_year/timevarying_pollution_severity_results.csv", c("organ", "outcome", "pollutant")),
  aj = list("prior_year_pollution_quartile_aalen_johansen_cif/prior_year_pollution_quartile_cif_at_1_3_5_10_years.csv", c("organ", "pollutant", "quartile", "time"))
)
summaries <- list()
for (stem in c("heart_durable_lvad", "kidney_dialysis_vintage")) {
  filename <- paste0(stem, "_aalen_johansen_cif_cif_at_1_3_5_10_years.csv")
  specs[[paste0("aj_", stem)]] <- list(
    file.path("figures/organ_subgroup_aalen_johansen_cif", filename),
    c("pollutant", "subgroup", "quartile", "time"),
    file.path("figures/organ_subgroup_aalen_johansen_cif/rolling_365d_20260922", filename)
  )
}
for (name in names(specs)) {
  spec <- specs[[name]]
  before <- read_csv(file.path(run, "before", spec[[1]]), show_col_types = FALSE)
  after_path <- if (length(spec) >= 3L) spec[[3]] else spec[[1]]
  after <- read_csv(file.path("output", after_path), show_col_types = FALSE)
  keys <- spec[[2]]
  stopifnot(all(keys %in% names(before)), all(keys %in% names(after)),
             !anyDuplicated(before[keys]), !anyDuplicated(after[keys]))
  metrics <- intersect(c("n", "people", "candidate_organ_episodes", "intervals", "adverse_events",
                         "hazard_ratio", "estimate", "coefficient", "std_error", "hr_conf_low", "hr_conf_high", "conf_low", "conf_high", "p_value",
                         "cif_adverse", "cif_transplant_or_improvement", "n_risk"),
                       intersect(names(before), names(after)))
  comparison <- full_join(select(before, all_of(c(keys, metrics))),
                           select(after, all_of(c(keys, metrics))), by = keys,
                           suffix = c("_before", "_after"))
  for (metric in metrics) {
    comparison[[paste0(metric, "_change")]] <- comparison[[paste0(metric, "_after")]] - comparison[[paste0(metric, "_before")]]
  }
  if ("hazard_ratio" %in% metrics) {
    comparison <- comparison %>% mutate(
      hr_relative_change_percent = 100 * (hazard_ratio_after / hazard_ratio_before - 1),
      direction_changed = (hazard_ratio_before > 1) != (hazard_ratio_after > 1),
      significance_changed = (p_value_before < .05) != (p_value_after < .05)
    )
    summaries[[name]] <- tibble(analysis = name, rows = nrow(comparison),
                                 direction_changes = sum(comparison$direction_changed, na.rm = TRUE),
                                 significance_changes = sum(comparison$significance_changed, na.rm = TRUE),
                                 max_absolute_relative_hr_change_percent = max(abs(comparison$hr_relative_change_percent), na.rm = TRUE))
  }
  write_csv(comparison, file.path(run, paste0("comparison_", name, ".csv")))
}
write_csv(bind_rows(summaries), file.path(run, "comparison_summary.csv"))
print(bind_rows(summaries))

table_path <- "tables/table1_baseline_waitlist_characteristics/table1_baseline_characteristics_by_organ.csv"
before <- read_csv(file.path(run, "before", table_path), show_col_types = FALSE)
after <- read_csv(file.path("output", table_path), show_col_types = FALSE)
writeLines(trimws(capture.output(all.equal(before, after)), which = "right"), file.path(run, "table1_comparison.txt"))
non_exposure <- function(x) x[!grepl("pollution exposure", x[[1]], ignore.case = TRUE), ]
stopifnot(identical(non_exposure(before), non_exposure(after)))

report <- c("# Rolling exposure revision results", "",
  "Baseline PM2.5 and O3 now use the 365 days before listing. NO2 uses the 12 complete calendar months before the listing month, day weighted, with annual approximations before 2019. Daily PM2.5 coverage also extends from 2021 to 2024. These comparisons therefore combine changes in exposure definition and eligibility; they do not isolate one change at a time.", "")
for (analysis in c("baseline", "timevarying")) {
  dat <- read_csv(file.path(run, paste0("comparison_", analysis, ".csv")), show_col_types = FALSE)
  pollutant <- if (analysis == "baseline") dat$exposure else dat$pollutant
  report <- c(report, paste0("## ", analysis, " primary Cox models"), "",
                "| Organ | Pollutant | Previous HR | Updated HR (95% CI) | Previous N | Updated N |",
                "|---|---|---:|---:|---:|---:|")
  n_col <- if (analysis == "baseline") "n" else "candidate_organ_episodes"
  for (i in seq_len(nrow(dat))) {
    report <- c(report, sprintf("| %s | %s | %.3f | %.3f (%.3f-%.3f) | %s | %s |",
                                dat$organ[i], pollutant[i], dat$hazard_ratio_before[i],
                                dat$hazard_ratio_after[i], dat$conf_low_after[i], dat$conf_high_after[i],
                                format(dat[[paste0(n_col, "_before")]][i], big.mark = ","),
                                format(dat[[paste0(n_col, "_after")]][i], big.mark = ",")))
  }
  report <- c(report, "")
}
aj <- read_csv(file.path(run, "comparison_aj.csv"), show_col_types = FALSE) %>%
  filter(pollutant %in% c("PM2.5", "NO2"), quartile %in% c("Q1 lowest", "Q4 highest"),
         (organ == "Heart" & time == 1) | (organ == "Lung" & time == 3) |
           (organ == "Liver" & time == 3) | (organ == "Kidney" & time == 10))
report <- c(report, "## AJ cumulative incidence of death/deterioration", "",
              "| Organ | Years | Pollutant | Quartile | Previous (%) | Updated (%) |",
              "|---|---:|---|---|---:|---:|")
for (i in seq_len(nrow(aj))) {
  report <- c(report, sprintf("| %s | %s | %s | %s | %.2f | %.2f |",
                              aj$organ[i], aj$time[i], aj$pollutant[i], aj$quartile[i],
                              100 * aj$cif_adverse_before[i], 100 * aj$cif_adverse_after[i]))
}
sens <- read_csv(file.path(run, "comparison_sensitivities.csv"), show_col_types = FALSE) %>%
  filter(analysis == "Baseline cause-specific Cox", organ %in% c("HR", "LU"),
         model == "Multipollutant model: PM2.5 + NO2")
report <- c(report, "", "## Baseline two-pollutant thoracic models", "",
              "| Organ | Term | Previous HR | Updated HR (95% CI) | Updated P |",
              "|---|---|---:|---:|---:|")
for (i in seq_len(nrow(sens))) {
  report <- c(report, sprintf("| %s | %s | %.3f | %.3f (%.3f-%.3f) | %.4g |",
                              sens$organ[i], sens$exposure[i], sens$hazard_ratio_before[i],
                              sens$hazard_ratio_after[i], sens$conf_low_after[i],
                              sens$conf_high_after[i], sens$p_value_after[i]))
}
report <- c(report, "", "## Scope and interpretation", "")
report <- c(report, "PM2.5 HRs are per 5 ug/m3; NO2 and O3 HRs are per 10 ppb.", "",
              "Scope: current baseline and monthly time-varying models, all three baseline adjustment tiers, subgroup and sensitivity analyses, severity associations, AJ curves, and their current figures/tables were regenerated. Superseded whole-waitlist-average, KM, and Fine-Gray analyses remain archived rather than being relabeled as current results.", "",
              "Table 1: candidate counts, demographics, clinical characteristics, and outcome rows are unchanged. The three exposure rows now describe prelisting exposures rather than whole-waitlist averages.", "",
              "The very large relative changes in the full-coefficient audit concern the kidney missing-race indicator, whose previous HR was near zero; these are nuisance-category coefficients, not pollutant associations. All pollutant estimates are separately compared in the baseline, model-tier, time-varying, subgroup, and sensitivity files.", "",
              "Figure QA: 48 updated PNG exports passed image checks and 48 updated PDFs opened successfully. Contact sheets are in figure_qa/. The main forest PDF was rendered separately for visual inspection. The full AJ and LVAD/dialysis AJ sets are in rolling_365d_20260922 subfolders to avoid locked previous exports.", "",
              "Interpretation caution: the existing map calculation uses 1 - exp(-H_death(365)) from a cause-specific Cox model. This is a net-risk transformation, not the cumulative incidence of death/deterioration before transplant. The latter requires the competing-event hazard as well. This exposure revision preserves the map calculation; it does not resolve that separate estimand issue.", "",
              "All model tiers, full coefficients, subgroups, sensitivities, severity associations, and AJ horizons have separate comparison CSVs in this directory. A changed significance classification alone does not establish a statistically significant difference between the old and new estimates.")
writeLines(report, file.path(run, "results_comparison.md"))
