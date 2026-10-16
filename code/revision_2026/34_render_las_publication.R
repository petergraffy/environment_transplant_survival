source("code/revision_2026/common.R")
dest <- file.path(revision_dir,"las_era_lung","render_logs")
dir.create(dest, recursive=TRUE,showWarnings=FALSE)
scripts <- c(
  "code/figures_tables/82_plot_figure1_prior_year_and_timevarying_cox.R",
  "code/figures_tables/86_plot_prior_year_and_timevarying_pollution_subgroup_forest.R",
  "code/figures_tables/87_plot_heart_lung_pm25_no2_subgroup_forest.R",
  "code/figures_tables/83_plot_prior_year_aalen_johansen_simplified.R",
  "code/revision_2026/05_aj_risk_tables.R",
  "code/figures_tables/69_make_table1_baseline_characteristics.R",
  "code/figures_tables/75_make_organ_score_formula_table.R",
  "code/figures_tables/89_make_updated_sensitivity_wide_tables.R",
  "code/figures_tables/92_plot_pollution_absolute_1yr_survival_maps.R",
  "code/figures_tables/93_plot_candidates_heart_lung_absolute_risk_abc.R",
  "code/figures_tables/94_plot_transplant_centers_by_organ_four_panel.R",
  "code/figures_tables/61_plot_waitlist_candidates_zcta_map.R",
  "code/figures_tables/62_plot_waitlist_geography_ab_maps.R")
args <- commandArgs(TRUE)
if(length(args)) scripts <- scripts[as.integer(args)]
for(script in scripts) {
  log_revision("Rendering ",script)
  path <- file.path(dest,paste0(basename(script),".log"))
  status <- system2(file.path(R.home("bin"),"Rscript.exe"),shQuote(script),stdout=path,stderr=path)
  if(status != 0) {cat(readLines(path),sep="\n"); stop("Render failed: ",script)}
  log_revision("Completed ",script)
}
