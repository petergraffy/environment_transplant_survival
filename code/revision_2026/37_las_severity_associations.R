source("code/revision_2026/common.R")
library(fixest)
setDTthreads(2)
b <- as_tibble(readRDS(file.path(revision_dir,"cache","las_era_lung","prepared.rds"))$cohort)
d <- readRDS(file.path(revision_dir,"cache","las_era_lung","monthly_intervals.rds"))
e <- load_script_prefix("code/80_prior_year_pollution_baseline_severity_associations.R","analysis_dat")
labels <- c(pm25="PM2.5 per 5 ug/m3",no2="NO2 per 10 ppb",o3="O3 per 10 ppb")
rows <- list(); coefs <- list()
for(p in names(pollutants)) {
  rows[[p]] <- e$fit_severity_model(b,"LU","baseline_las","Official LAS at listing",p,pollutants[[p]],labels[[p]],"gaussian")
  coefs[[p]] <- e$fit_severity_model_coefficients(b,"LU","baseline_las","Official LAS at listing",p,pollutants[[p]],labels[[p]],"gaussian")
}
update <- function(path,rows,key="organ",lung="LU") {
  old <- fread(path)
  out <- rbindlist(list(old[get(key)!=lung],as.data.table(rows)),fill=TRUE)
  fwrite(out,path); out
}
dir <- "output/prior_year_pollution_baseline_severity_associations"
res <- update(file.path(dir,"prior_year_pollution_baseline_severity_associations.csv"),bind_rows(rows))
fwrite(res[,estimate_ci:=sprintf("%.2f (%.2f-%.2f)",estimate,conf_low,conf_high)],file.path(dir,"prior_year_pollution_baseline_severity_associations_table.csv"))
full <- update(file.path(dir,"prior_year_pollution_baseline_severity_full_coefficients_all_terms.csv"),bind_rows(coefs))
fwrite(full[term_type!="listing_center_fixed_effect"],file.path(dir,"prior_year_pollution_baseline_severity_full_coefficients_supplement.csv"))
e <- load_script_prefix("code/84_timevarying_pollution_severity_models.R","analysis_dat")
e$severity_outcome <- function(org) list(term="las_observed",label="Official LAS (points)")
rows <- list()
for(p in names(pollutants)) {
  term <- c(pm25="pm25_interval_5ug",no2="no2_interval_10ppb",o3="o3_interval_10ppb")[[p]]
  rows[[p]] <- e$fit_tv_severity_model(d,"LU",p,term,labels[[p]])
}
dir <- "output/timevarying_pollution_severity_models_with_listing_year"
res <- update(file.path(dir,"timevarying_pollution_severity_results.csv"),bind_rows(rows))
fwrite(res[,estimate_ci:=sprintf("%.2f (%.2f-%.2f)",estimate,conf_low,conf_high)],file.path(dir,"timevarying_pollution_severity_table.csv"))
writeLines(c("Lung severity response is observed official LAS in points; missing outcomes are not imputed.",
  "Lung listing and follow-up end March 8, 2023. Baseline response models omit listing year; time-varying response models include listing year, preserving their existing adjustment sets.",
  "Time-varying response regressions retain the existing interval-level fixed-effect least-squares specification and its model-based uncertainty; they are descriptive associations, not causal mediation estimates."),
  file.path(revision_dir,"las_era_lung","severity_response_notes.txt"))
