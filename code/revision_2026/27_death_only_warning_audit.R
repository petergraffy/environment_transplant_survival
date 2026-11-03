source("code/revision_2026/common.R")
source("code/revision_2026/outcome_components.R")
setDTthreads(2)
dest <- file.path(revision_dir, "death_only_sensitivity")
r <- fread(file.path(dest, "pollutant_results.csv"))
b <- baseline_data()
components <- primary_outcome_components(b)
b <- merge(b[WL_ORG == "LU"], components[, .(PERS_ID, WL_ORG, death_event)], by = c("PERS_ID", "WL_ORG"))
audits <- list()
audit_fit <- function(d, exposure_term, rhs, model_type, p, tv = FALSE) {
  ix <- r$model == model_type & r$organ == "LU" & r$pollutant == p & r$endpoint == "Death only"
  original <- r[ix]
  stopifnot(nrow(original) == 1, grepl("variable  6", original$warnings))
  f <- as.formula(paste(if(tv) "Surv(tstart,tstop,analysis_event)~" else "Surv(followup_days,analysis_event)~",
    paste(c(exposure_term, rhs), collapse = "+")))
  d[, analysis_event := death_event * if(tv) tv_adverse_event else adverse_event]
  missing <- d[as.character(race) == "Missing"]
  x <- droplevels(as.data.frame(d[as.character(race) != "Missing"]))
  warnings <- character()
  fit <- withCallingHandlers(coxph(f, data = x, ties = "efron", robust = FALSE, x = FALSE, y = FALSE),
    warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")})
  z <- as.data.table(broom::tidy(fit, conf.int = TRUE))[term == exposure_term]
  audits[[length(audits) + 1L]] <<- data.table(model = model_type, organ = "LU", pollutant = p,
    warning_term = "raceMissing", missing_race_candidates = uniqueN(missing$waitlist_row_id),
    missing_race_intervals = nrow(missing), missing_race_deaths = sum(missing$analysis_event),
    full_cohort_HR = original$hazard_ratio, diagnostic_HR = exp(z$estimate),
    full_cohort_SE = original$std.error, diagnostic_SE = z$std.error,
    HR_absolute_difference = abs(original$hazard_ratio - exp(z$estimate)),
    diagnostic_warnings = paste(warnings, collapse = "; "))
  stopifnot(sum(missing$analysis_event) == 0, length(warnings) == 0,
    abs(original$hazard_ratio - exp(z$estimate)) < 1e-4)
  fwrite(rbindlist(audits), file.path(dest, "lung_warning_audit.csv"))
}
for (p in c("pm25", "o3")) audit_fit(as.data.table(eligible(b, "LU", pollutants[[p]])),
  pollutants[[p]], baseline_adjustment, "Baseline", p)
env <- load_script_prefix("code/79_timevarying_pollution_cox_svi.R", "analysis_dat")
ex <- env$read_exposure_tables()
for (p in c("pm25", "o3")) {
  d <- env$make_month_intervals(b, daily_pollution_end_date)
  d <- env$add_time_updated_scores(d, list(liver = data.table(), heart = data.table()))
  d <- ex[[p]][d, on = c("zip" = "candidate_zip", "year", "month")]
  term <- if(p == "pm25") "pm25_interval_5ug" else "o3_interval_10ppb"
  d[, (term) := if(p == "pm25") pm25_interval_ug_m3 / 5 else o3_interval_ppb / 10]
  d[b[, .(waitlist_row_id, death_event)], on = "waitlist_row_id", death_event := i.death_event]
  adj <- env$primary_terms("LU")
  needed <- c("tstart", "tstop", "tv_adverse_event", term, adj$vars)
  d <- d[complete.cases(d[, ..needed]) & tstop > tstart]
  audit_fit(d, term, adj$terms, "Time-varying", p, TRUE)
}
print(rbindlist(audits))
