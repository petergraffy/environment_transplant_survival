# Preserve the transition-censored bundle and fit the same models with longer follow-up.
options(las.bundle = "las_extended_followup", las.extend_followup = TRUE)
source("code/revision_2026/30_prepare_las_era_lung.R")
source("code/revision_2026/31_fit_las_era_lung.R")
transition <- as.Date("2023-03-09")
audit <- rbindlist(lapply(c("pm25", "no2", "o3"), function(p) {
  term <- terms[[p]]
  end <- if(p == "no2") as.Date("2025-12-31") else as.Date("2024-12-31")
  x <- d[interval_end_date <= end & complete.cases(d[, c(term, "age", "sex", "race", "zcta_svi_proxy", "listing_center"), with = FALSE])]
  post_days <- pmax(0, as.numeric(x$interval_end_date - pmax(x$interval_start_date, transition)) + 1)
  data.table(pollutant = p, N = uniqueN(x$PERS_ID), person_days = sum(x$tstop-x$tstart),
    post_transition_people = uniqueN(x[interval_end_date >= transition, PERS_ID]),
    post_transition_days = sum(post_days), post_transition_percent = 100*sum(post_days)/sum(x$tstop-x$tstart),
    adverse_events = sum(x$tv_adverse_event),
    post_transition_adverse_events = sum(x$tv_adverse_event[x$interval_end_date >= transition]),
    post_transition_days_with_observed_las = sum(post_days[x$las_unavailable == 0]))
}))
fwrite(audit, file.path(dest, "post_transition_followup.csv"))
old <- fread(file.path(revision_dir, "las_era_lung", "pollutant_results.csv"))
old[, followup := "Transition censored"]
new <- fread(file.path(dest, "pollutant_results.csv"))
new[, followup := "Extended; last LAS carried forward"]
comparison <- rbindlist(list(old, new), fill = TRUE)[variant == "Primary"]
fwrite(comparison, file.path(dest, "primary_comparison.csv"))
stopifnot(nrow(b) == 42874L, all(prep$scores$score_date <= as.Date("2023-03-08")),
  all(d$las_source_date <= d$start_date, na.rm = TRUE),
  all(d[start_date >= transition, las_source_date] < transition, na.rm = TRUE))
writeLines(c("Enrollment ends March 8, 2023. Outcomes followed through December 31, 2025, subject to observed follow-up.",
  "LAS updates use only measurements available by March 8, 2023; the last observed score is carried forward without future backfill.",
  "Baseline prior-year exposure remains fixed; time-varying PM2.5/O3 follow-up ends December 31, 2024, NO2 December 31, 2025.",
  "Multipollutant time-varying models use the common December 31, 2024 endpoint.",
  "Missing and pediatric LAS handling is unchanged. Post-transition LAS is last documented severity, not updated severity.",
  "The transition-censored results and publication figures remain preserved in their existing locations."),
  file.path(dest, "methods.txt"))
print(comparison[, .(followup, model, term, N, events, HR, HR_low, HR_high)])
print(audit)
