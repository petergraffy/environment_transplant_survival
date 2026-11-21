#!/usr/bin/env Rscript

# Run from the repository root. Each script runs in a fresh R process.
steps <- list(
  inputs = c(
    "code/prepare_publication_inputs.R"
  ),
  cohort = c(
    "code/50_primary_waitlist_period_pollution_cox.R",
    "code/prepare_daily_pollution_inputs.R",
    "code/prepare_rolling_prior_pollution.R",
    "code/revision_2026/02_prepare_registration_sensitivities.R",
    "code/figures_tables/19_map_transplant_centers.R"
  ),
  models = c(
    "code/81_prior_year_pollution_cox_svi.R",
    "code/79_timevarying_pollution_cox_svi.R",
    "code/80_prior_year_pollution_baseline_severity_associations.R",
    "code/84_timevarying_pollution_severity_models.R",
    "code/78_prior_year_pollution_quartile_aalen_johansen_cif.R",
    "code/86_prior_year_and_timevarying_pollution_subgroup_cox.R",
    "code/88_updated_methods_sensitivity_models.R"
  ),
  sensitivities = c(
    "code/revision_2026/03_baseline_diagnostics_eras.R",
    "code/revision_2026/04_svi_validation.R",
    "code/revision_2026/07_registration_models.R",
    "code/revision_2026/08_timevarying_diagnostics.R",
    "code/revision_2026/13_thoracic_warning_audit.R",
    "code/revision_2026/14_baseline_warning_audit.R",
    "code/revision_2026/15_kidney_late_band_audit.R",
    "code/revision_2026/16_pre_covid_competing_risks.R",
    "code/revision_2026/25_death_only_sensitivity.R",
    "code/revision_2026/27_death_only_warning_audit.R",
    "code/revision_2026/28_retransplant_candidate_counts.R"
  ),
  publication = c(
    "code/revision_2026/29_audit_official_lung_scores.R",
    "code/revision_2026/30_prepare_las_era_lung.R",
    "code/revision_2026/31_fit_las_era_lung.R",
    "code/revision_2026/32_sync_las_publication_inputs.R",
    "code/revision_2026/35_las_lung_diagnostics_sensitivities.R",
    "code/90_prior_year_initial_listing_cox_model_set.R",
    "code/revision_2026/37_las_severity_associations.R",
    "code/revision_2026/40_las_extended_followup.R",
    "code/revision_2026/43_transplant_improvement_baseline.R"
  ),
  displays = c(
    "code/figures_tables/82_plot_figure1_prior_year_and_timevarying_cox.R",
    "code/figures_tables/83_plot_prior_year_aalen_johansen_simplified.R",
    "code/revision_2026/05_aj_risk_tables.R",
    "code/figures_tables/71_plot_heart_lvad_kidney_dialysis_aalen_johansen_cif.R",
    "code/figures_tables/85_plot_pollution_severity_associations.R",
    "code/figures_tables/86_plot_prior_year_and_timevarying_pollution_subgroup_forest.R",
    "code/figures_tables/87_plot_heart_lung_pm25_no2_subgroup_forest.R",
    "code/figures_tables/91_plot_initial_listing_model_tier_forest.R",
    "code/figures_tables/60_plot_primary_pollution_study_period_maps.R",
    "code/figures_tables/92_plot_pollution_absolute_1yr_survival_maps.R",
    "code/figures_tables/93_plot_candidates_heart_lung_absolute_risk_abc.R",
    "code/figures_tables/94_plot_transplant_centers_by_organ_four_panel.R",
    "code/figures_tables/61_plot_waitlist_candidates_zcta_map.R",
    "code/figures_tables/62_plot_waitlist_geography_ab_maps.R",
    "code/revision_2026/42_cohort_flow_original_layout.R",
    "code/figures_tables/69_make_table1_baseline_characteristics.R",
    "code/figures_tables/75_make_organ_score_formula_table.R",
    "code/figures_tables/77_make_svi_proxy_supplement_table.R",
    "code/figures_tables/89_make_updated_sensitivity_wide_tables.R",
    "code/revision_2026/33_las_publication_tables_flow.R",
    "code/revision_2026/10_supplement_tables.R",
    "code/revision_2026/11_validation_and_diagnostic_figures.R",
    "code/revision_2026/17_svi_validation_figure.R",
    "code/revision_2026/26_death_only_report.R",
    "code/revision_2026/41_las_extended_report.R"
  ),
  validate = c("code/revision_2026/39_validate_las_publication.R"),
  qa = c("code/revision_2026/38_las_publication_visual_qa.py")
)

args <- commandArgs(trailingOnly = TRUE)
scripts <- unlist(steps, use.names = FALSE)
missing <- scripts[!file.exists(scripts)]
if (length(missing)) stop("Run from the repository root. Missing pipeline scripts: ", paste(missing, collapse = ", "))
if (identical(args, "--list")) {
  for (stage in names(steps)) cat(stage, ":\n", paste0("  ", steps[[stage]], collapse = "\n"), "\n", sep = "")
  quit(save = "no")
}
if (identical(args, "--check")) {
  invisible(lapply(list.files("code", pattern = "[.]R$", recursive = TRUE, full.names = TRUE), parse))
  message("Pipeline scripts exist and R sources parse.")
  quit(save = "no")
}
if (any(!args %in% c(names(steps), "--download"))) stop("Usage: Rscript code/run_paper.R [--list|--check|stage ...] [--download]")
stages <- setdiff(args, "--download")
if (!length(stages)) stages <- setdiff(names(steps), "qa")
stages <- intersect(names(steps), stages)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
python <- Sys.getenv("PAPER_PYTHON", Sys.which(if (.Platform$OS.type == "windows") "py" else "python3"))
log_dir <- file.path("output", "paper_run_logs", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
for (stage in stages) for (script in steps[[stage]]) {
  log <- file.path(log_dir, paste0(basename(script), ".log"))
  exe <- if (endsWith(script, ".py")) python else rscript
  extra <- if (stage == "inputs" && "--download" %in% args) "--download" else character()
  if (!nzchar(exe)) stop("Set PAPER_PYTHON to the Python executable for figure QA.")
  message(format(Sys.time(), "%H:%M:%S"), " [", stage, "] ", script)
  status <- system2(exe, shQuote(c(script, extra)), stdout = log, stderr = log)
  if (status != 0L) stop("Failed: ", script, ". See ", log, call. = FALSE)
}
message("Completed: ", paste(stages, collapse = ", "), ". Logs: ", log_dir)
