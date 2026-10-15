source("code/revision_2026/common.R")
setDTthreads(2)
dest <- file.path(revision_dir, "las_era_lung")
prep <- readRDS(file.path(revision_dir, "cache", "las_era_lung", "prepared.rds"))
b <- prep$cohort
b[, `:=`(age_group = cut(age, c(-Inf, 17, 39, 59, Inf), labels = c("<18", "18-39", "40-59", "60+")),
  race_group = fcase(as.character(race) == "WHITE", "White", as.character(race) == "BLACK", "Black",
    as.character(race) == "ASIAN", "Asian", default = "Other/Unknown"))]
results <- list(); coefficients <- list(); logs <- list()
fit_one <- function(d, terms, model_type, variant = "Primary", subgroup = "All",
                    adjust = c("age", "sex", "race", "zcta_svi_proxy"), center = TRUE,
                    severity = TRUE, death_only = FALSE) {
  tv <- model_type == "Time-varying"
  event <- if(tv) "tv_adverse_event" else "adverse_event"
  severity_terms <- if(tv && severity) c("las_per10", "las_unavailable", "las_under12") else character()
  if(variant == "Previous proxy") severity_terms <- "organ_score_tv"
  needed <- unique(c(terms, adjust, severity_terms, "listing_center", event,
    if(tv) c("tstart", "tstop") else "followup_days"))
  x <- copy(d[complete.cases(d[, ..needed])])
  x <- if(tv) x[tstop > tstart] else x[followup_days > 0]
  x[, analysis_event := as.integer(get(event) == 1 & (!death_only | death_event == 1))]
  rhs <- c(adjust, severity_terms)
  rhs <- rhs[vapply(rhs, function(v) uniqueN(x[[v]]) > 1, logical(1))]
  f <- as.formula(paste(if(tv) "Surv(tstart,tstop,analysis_event)~" else "Surv(followup_days,analysis_event)~",
    paste(c(terms, rhs, if(center) "strata(listing_center)"), collapse = "+")))
  log_revision(model_type, " / ", variant, " / ", subgroup, " / ", paste(terms, collapse = "+"))
  warnings <- character()
  fit <- if(sum(x$analysis_event) < 50) simpleError("Fewer than 50 events") else
    tryCatch(withCallingHandlers(coxph(f, data = droplevels(as.data.frame(x)), ties = "efron", robust = FALSE,
      x = FALSE, y = FALSE), warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}), error = function(e) e)
  status <- if(inherits(fit, "error")) conditionMessage(fit) else if(length(warnings)) paste(unique(warnings), collapse = "; ") else "OK"
  meta <- list(model = model_type, variant = variant, subgroup = subgroup,
    N = uniqueN(x$waitlist_row_id), intervals = nrow(x), events = sum(x$analysis_event),
    formula = paste(deparse(f), collapse = " "), status = status)
  logs[[length(logs) + 1L]] <<- as.data.table(meta)
  if(!inherits(fit, "error")) {
    z <- as.data.table(broom::tidy(fit, conf.int = TRUE))
    for(nm in names(meta)) set(z, j = nm, value = meta[[nm]])
    z[, `:=`(HR = exp(estimate), HR_low = exp(conf.low), HR_high = exp(conf.high))]
    coefficients[[length(coefficients) + 1L]] <<- z
    results[[length(results) + 1L]] <<- z[term %in% terms]
  }
  fwrite(rbindlist(results), file.path(dest, "pollutant_results.csv"))
  fwrite(rbindlist(coefficients), file.path(dest, "full_coefficients.csv"))
  fwrite(rbindlist(logs), file.path(dest, "model_log.csv"))
}
fit_suite <- function(d, terms, model_type) {
  for(p in names(terms)) {
    term <- terms[[p]]
    fit_one(d, term, model_type)
    fit_one(d, term, model_type, "No center strata", center = FALSE)
    fit_one(d[multi_organ_candidate == 0], term, model_type, "Single organ only")
    fit_one(d, term, model_type, "Adjust multi-organ", adjust = c("age", "sex", "race", "zcta_svi_proxy", "multi_organ_candidate"))
    fit_one(d, term, model_type, "Death only", death_only = TRUE)
    if(model_type == "Time-varying") {
      fit_one(d, term, model_type, "Previous proxy")
      fit_one(d, term, model_type, "No severity adjustment", severity = FALSE)
      fit_one(d[las_unavailable == 0], term, model_type, "Observed LAS only")
    }
    for(g in levels(b$age_group)) fit_one(d[age_group == g], term, model_type, "Age subgroup", g,
      adjust = c("sex", "race_group", "zcta_svi_proxy"))
    for(g in c("F", "M")) fit_one(d[as.character(sex) == g], term, model_type, "Sex subgroup", if(g == "F") "Female" else "Male",
      adjust = c(if(model_type == "Time-varying") "age_interval" else "age", "race_group", "zcta_svi_proxy"))
    for(g in c("White", "Black", "Asian")) fit_one(d[race_group == g], term, model_type, "Race subgroup", g,
      adjust = c(if(model_type == "Time-varying") "age_interval" else "age", "sex", "zcta_svi_proxy"))
  }
  fit_one(d, terms[c("pm25", "no2")], model_type, "PM2.5 + NO2")
  fit_one(d, terms, model_type, "PM2.5 + NO2 + O3")
}
fit_suite(b, pollutants, "Baseline")
env <- load_script_prefix("code/79_timevarying_pollution_cox_svi.R", "analysis_dat")
d <- env$make_month_intervals(b, prep$cutoff)
d[, `:=`(start_date = index_date + tstart, age_interval = age + tstart / 365.25)]
d[b[, .(waitlist_row_id, death_event)], on = "waitlist_row_id", death_event := i.death_event]
setkey(prep$scores, PERS_ID, score_date)
d[, las_observed := prep$scores[.SD, on = .(PERS_ID, score_date = start_date), roll = Inf, x.las]]
d[, las_source_date := prep$scores[.SD, on = .(PERS_ID, score_date = start_date), roll = Inf, x.score_date]]
d[, las_under12 := as.integer(age_interval < 12)]
d[las_under12 == 1, las_observed := NA_real_]
d[, las_unavailable := as.integer(!is.finite(las_observed))]
d[, las_per10 := fifelse(las_unavailable == 1, prep$imputation, las_observed) / 10]
stopifnot(all(d$las_source_date <= d$start_date, na.rm = TRUE),
  all(d$interval_end_date <= prep$cutoff), sum(d$tv_adverse_event) == sum(b$adverse_event))
coverage <- d[, .(intervals = .N, people = uniqueN(PERS_ID), person_days = sum(tstop - tstart),
  observed_las_days = sum((tstop - tstart)[las_unavailable == 0]),
  zero_las_days = sum((tstop - tstart)[las_observed == 0], na.rm = TRUE),
  under12_days = sum((tstop - tstart)[las_under12 == 1]),
  missing_las_event_count = sum(tv_adverse_event[las_unavailable == 1]))]
fwrite(coverage, file.path(dest, "time_updated_las_coverage.csv"))
ex <- env$read_exposure_tables()
d <- ex$pm25[d, on = c("zip" = "candidate_zip", "year", "month")]
d <- ex$o3[d, on = c("zip", "year", "month")]
d <- ex$no2_monthly[d, on = c("zip", "year", "month")]
d <- ex$no2_annual[d, on = c("zip", "year")]
d[, `:=`(pm25_interval_5ug = pm25_interval_ug_m3 / 5, o3_interval_10ppb = o3_interval_ppb / 10,
  no2_interval_10ppb = fifelse(!is.na(no2_interval_ppb), no2_interval_ppb, no2_annual_ppb) / 10)]
terms <- c(pm25 = "pm25_interval_5ug", no2 = "no2_interval_10ppb", o3 = "o3_interval_10ppb")
saveRDS(d, file.path(revision_dir, "cache", "las_era_lung", "monthly_intervals.rds"), compress = FALSE)
fit_suite(d, terms, "Time-varying")
r <- rbindlist(results)
print(r[variant == "Primary", .(model, term, N, events, HR, HR_low, HR_high, status)])
log_revision("LAS-era lung fits complete")
