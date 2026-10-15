source("code/revision_2026/common.R")
source("code/revision_2026/outcome_components.R")
setDTthreads(2)
dest <- file.path(revision_dir, "death_only_sensitivity")
dir.create(dest, showWarnings = FALSE)
b <- baseline_data()
components <- primary_outcome_components(b)
b <- merge(b, components[, .(PERS_ID, WL_ORG, component, death_event)], by = c("PERS_ID", "WL_ORG"))
counts <- components[, .(N = .N, death = sum(component == "death"),
  deterioration = sum(component == "deterioration"), transplant = sum(component == "transplant"),
  improvement = sum(component == "improvement"), other_exit = sum(component == "other_exit"),
  administrative_censor = sum(component == "administrative_censor"),
  transplant_procedure_deaths = sum(terminal_removal_code %in% c(21, 23) & component == "transplant")), by = WL_ORG]
stopifnot(all(counts$N == rowSums(counts[, .(death, deterioration, transplant, improvement, other_exit, administrative_censor)])))
fwrite(counts, file.path(dest, "outcome_components_by_organ.csv"))
fwrite(components[, .N, by = .(WL_ORG, component, terminal_removal_code)], file.path(dest, "removal_code_audit.csv"))
primary <- rbindlist(list(
  fread("output/prior_year_pollution_cox_svi/prior_year_pollution_cox_svi_results.csv")[, `:=`(model = "Baseline", pollutant = exposure, candidate_organ_episodes = n)],
  fread("output/timevarying_pollution_cox_svi/timevarying_pollution_cox_svi_results.csv")[, model := "Time-varying"]), fill = TRUE)

fit_pair <- function(d, org, p, exposure_term, adjustment, model) {
  path <- file.path(dest, paste0(tolower(model), "_", org, "_", p, "_coefficients.csv"))
  if (file.exists(path)) return(invisible(NULL))
  tv <- model == "Time-varying"
  event <- if (tv) "tv_adverse_event" else "adverse_event"
  d[, analysis_event := get(event)]
  f <- as.formula(paste(if (tv) "Surv(tstart,tstop,analysis_event) ~" else "Surv(followup_days,analysis_event) ~",
    paste(c(exposure_term, adjustment), collapse = " + ")))
  output <- list()
  for (endpoint in c("Composite", "Death only")) {
    d[, analysis_event := as.integer(get(event) == 1 & (endpoint == "Composite" | death_event == 1))]
    log_revision(model, " / ", org, " / ", p, " / ", endpoint, " rows=", nrow(d), " events=", sum(d$analysis_event))
    warnings <- character()
    fit <- withCallingHandlers(coxph(f, data = d, ties = "efron", robust = FALSE, x = FALSE, y = FALSE),
      warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")})
    z <- as.data.table(broom::tidy(fit, conf.int = TRUE))
    z[, `:=`(model = model, organ = org, pollutant = p, endpoint = endpoint,
      candidate_organ_episodes = uniqueN(d$waitlist_row_id), intervals = nrow(d), events = sum(d$analysis_event),
      hazard_ratio = exp(estimate), HR_low = exp(conf.low), HR_high = exp(conf.high),
      formula = paste(deparse(f), collapse = " "), warnings = paste(unique(warnings), collapse = "; "))]
    if (endpoint == "Composite") {
      ix <- primary$model == model & primary$organ == org & primary$pollutant == p
      ref <- primary[ix]
      current <- z[term == exposure_term]
      stopifnot(nrow(ref) == 1, nrow(current) == 1,
        current$candidate_organ_episodes == ref$candidate_organ_episodes,
        current$events == ref$adverse_events,
        abs(current$hazard_ratio - ref$hazard_ratio) < 1e-5)
    }
    output[[endpoint]] <- z
    rm(fit); gc()
  }
  fwrite(rbindlist(output), path)
}

for (org in names(organ_names)) for (p in names(pollutants)) {
  d <- as.data.table(eligible(b, org, pollutants[[p]]))
  fit_pair(d, org, p, pollutants[[p]], baseline_adjustment, "Baseline")
}
rm(d); gc()
env <- load_script_prefix("code/79_timevarying_pollution_cox_svi.R", "analysis_dat")
ex <- env$read_exposure_tables()
scores <- env$make_time_updated_score_sources()
terms <- c(pm25 = "pm25_interval_5ug", no2 = "no2_interval_10ppb", o3 = "o3_interval_10ppb")
for (org in names(organ_names)) for (p in names(terms)) {
  path <- file.path(dest, paste0("time-varying_", org, "_", p, "_coefficients.csv"))
  if (file.exists(path)) next
  log_revision("Preparing monthly intervals: ", org, " / ", p)
  d <- env$make_month_intervals(b[WL_ORG == org], if(p == "no2") as.Date("2025-12-31") else daily_pollution_end_date)
  d <- env$add_time_updated_scores(d, scores)
  if (p == "pm25") {
    d <- ex$pm25[d, on = c("zip" = "candidate_zip", "year", "month")]
    d[, pm25_interval_5ug := pm25_interval_ug_m3 / 5]
  } else if (p == "o3") {
    d <- ex$o3[d, on = c("zip" = "candidate_zip", "year", "month")]
    d[, o3_interval_10ppb := o3_interval_ppb / 10]
  } else {
    d <- ex$no2_monthly[d, on = c("zip" = "candidate_zip", "year", "month")]
    d <- ex$no2_annual[d, on = c("zip", "year")]
    d[, no2_interval_10ppb := fifelse(!is.na(no2_interval_ppb), no2_interval_ppb, no2_annual_ppb) / 10]
  }
  d[b[, .(waitlist_row_id, death_event)], on = "waitlist_row_id", death_event := i.death_event]
  adj <- env$primary_terms(org)
  needed <- c("tstart", "tstop", "tv_adverse_event", terms[[p]], adj$vars)
  keep <- c(needed, "waitlist_row_id", "death_event")
  d <- d[complete.cases(d[, ..needed]) & tstop > tstart, ..keep]
  gc()
  fit_pair(d, org, p, terms[[p]], adj$terms, "Time-varying")
  rm(d); gc()
}
all <- rbindlist(lapply(list.files(dest, pattern = "_coefficients.csv$", full.names = TRUE), fread))
fwrite(all, file.path(dest, "full_model_results.csv"))
results <- all[term %in% c(unname(pollutants), unname(terms))]
stopifnot(nrow(results) == 48, uniqueN(results[, .(model, organ, pollutant, endpoint)]) == 48)
fwrite(results, file.path(dest, "pollutant_results.csv"))
log_revision("Death-only sensitivity complete; all 24 composite refits matched saved primary results")
