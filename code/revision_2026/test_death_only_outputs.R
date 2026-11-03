source("code/revision_2026/common.R")
dest <- file.path(revision_dir, "death_only_sensitivity")
r <- fread(file.path(dest, "pollutant_results.csv"))
r[, warnings := fifelse(is.na(warnings), "", as.character(warnings))]
c <- fread(file.path(dest, "outcome_components_by_organ.csv"))
stopifnot(nrow(r) == 48, uniqueN(r[, .(model, organ, pollutant, endpoint)]) == 48,
  sum(c$N) == 1017540,
  all(rowSums(c[, .(death, deterioration, transplant, improvement, other_exit, administrative_censor)]) == c$N),
  all(abs(exp(r$estimate) - r$hazard_ratio) < 1e-10),
  all(abs(exp(r$estimate - qnorm(.975) * r$std.error) - r$HR_low) < 1e-8),
  all(abs(exp(r$estimate + qnorm(.975) * r$std.error) - r$HR_high) < 1e-8),
  !anyNA(r$hazard_ratio))
flagged <- r[nzchar(warnings)]
audit <- fread(file.path(dest, "lung_warning_audit.csv"))
stopifnot(nrow(flagged) == 4, all(flagged$organ == "LU"),
  all(flagged$endpoint == "Death only"), all(flagged$pollutant %in% c("pm25", "o3")),
  all(grepl("variable  6", flagged$warnings)), nrow(audit) == 4,
  all(audit$missing_race_deaths == 0), all(audit$HR_absolute_difference < 1e-4))
pairs <- merge(r[endpoint == "Composite"], r[endpoint == "Death only"],
  by = c("model", "organ", "pollutant"), suffixes = c("_primary", "_death"))
stopifnot(all(pairs$candidate_organ_episodes_primary == pairs$candidate_organ_episodes_death),
  all(pairs$intervals_primary == pairs$intervals_death), all(pairs$events_death < pairs$events_primary))
deaths <- merge(r[endpoint == "Death only"], c, by.x = "organ", by.y = "WL_ORG")
stopifnot(all(deaths$events <= deaths$death), all(deaths$candidate_organ_episodes <= deaths$N))
saved <- rbindlist(list(
  fread("output/prior_year_pollution_cox_svi/prior_year_pollution_cox_svi_results.csv")[,
    .(model = "Baseline", organ, pollutant = exposure, saved_N = n,
      saved_events = adverse_events, saved_HR = hazard_ratio, saved_low = conf_low, saved_high = conf_high)],
  fread("output/timevarying_pollution_cox_svi/timevarying_pollution_cox_svi_results.csv")[,
    .(model = "Time-varying", organ, pollutant, saved_N = candidate_organ_episodes,
      saved_events = adverse_events, saved_HR = hazard_ratio, saved_low = conf_low, saved_high = conf_high)]))
reproduced <- merge(r[endpoint == "Composite"], saved, by = c("model", "organ", "pollutant"))
stopifnot(nrow(reproduced) == 24, all(reproduced$candidate_organ_episodes == reproduced$saved_N),
  all(reproduced$events == reproduced$saved_events),
  all(abs(reproduced$hazard_ratio - reproduced$saved_HR) < 1e-5),
  all(abs(reproduced$HR_low - reproduced$saved_low) < 1e-5),
  all(abs(reproduced$HR_high - reproduced$saved_high) < 1e-5))
path <- "output/tables/table1_baseline_waitlist_characteristics/table1_baseline_characteristics_by_organ.csv"
old_text <- system2("git", c("show", paste0("HEAD:", path)), stdout = TRUE)
old <- fread(text = paste(old_text, collapse = "\n"))
new <- fread(path)
stopifnot(identical(old[section != "Waitlist outcome, n (%)"], new[section != "Waitlist outcome, n (%)"]))
writeLines(c("PASS: 24 paired sensitivities, 48 exposure estimates.",
  "PASS: all 24 primary composite HRs, CIs, candidate counts, and event counts reproduced.",
  "PASS: identical candidates and intervals within each model pair; death events are a strict subset.",
  "PASS: coefficients, standard errors, HRs and CIs reconcile.",
  "PASS: four lung warnings localized to zero-death missing-race category; diagnostic HR changes <0.0001. Full-cohort fits retained and flagged.",
  "PASS: six mutually exclusive outcome categories sum to each organ denominator and 1,017,540 total.",
  "PASS: all non-outcome Table 1 cells identical to HEAD."), file.path(dest, "validation.txt"))
message("All death-only output checks passed")
