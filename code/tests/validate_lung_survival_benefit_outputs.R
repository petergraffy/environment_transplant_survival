# Validate completed real-data outputs. Writes aggregate checks only.
source("code/r_runtime.R"); ensure_user_library()
suppressPackageStartupMessages({library(data.table); library(survival)})
source("code/lung_survival_benefit_helpers.R")
dest <- Sys.getenv("LUNG_BENEFIT_OUTPUT", "output/lung_pollution_survival_benefit")
prepared <- readRDS(file.path(dest, "cache", "prepared.rds"))
fitted <- readRDS(file.path(dest, "cache", "fitted.rds"))
checks <- list(); missing <- list()
for (phase in c("pre", "post")) {
  d <- prepared[[phase]]
  stopifnot(all(d$tstop > d$tstart), !anyDuplicated(d[, .(PERS_ID, tstart)]),
    all(d[, sum(event), by = PERS_ID]$V1 <= 1L),
    all(d$exposure_source_end < d$start_date, na.rm = TRUE),
    all(d$las_source_date <= d$start_date, na.rm = TRUE))
  if (phase == "post") stopifnot(all(d$las_source_date <= d$index_date, na.rm = TRUE))
  contiguous <- d[order(PERS_ID, tstart), .(ok = all(tstart[-1] == head(tstop, -1))), by = PERS_ID]
  stopifnot(all(contiguous$ok))
  checks[[phase]] <- data.table(phase, people = uniqueN(d$PERS_ID), intervals = nrow(d),
    deaths = sum(d$event), future_exposure_intervals = sum(d$exposure_source_end >= d$start_date, na.rm = TRUE),
    future_las_intervals = sum(d$las_source_date > d$start_date, na.rm = TRUE))
  for (p in c("pm25", "no2", "o3")) {
    m <- fitted$models[[p]][[paste0(phase, "transplant")]]
    stopifnot(sum(m$data$event) == m$fit$nevent, all(is.finite(coef(m$fit))))
    missing[[paste(phase, p)]] <- data.table(phase, pollutant = p,
      available_people = uniqueN(d$PERS_ID), modeled_people = uniqueN(m$data$PERS_ID),
      available_intervals = nrow(d), modeled_intervals = nrow(m$data),
      available_deaths = sum(d$event), modeled_deaths = sum(m$data$event),
      available_person_years = sum(d$tstop - d$tstart) / 365.25,
      modeled_person_years = sum(m$data$tstop - m$data$tstart) / 365.25)
    r <- fitted$reference[[p]]
    stopifnot(!anyDuplicated(r$PERS_ID),
      nrow(lung_supported_reference(r, fitted$models[[p]]$pretransplant,
        fitted$models[[p]]$posttransplant, unique(fitted$point$horizon_days))) == nrow(r))
  }
}
c <- fread(file.path(dest, "survival_benefit_curves.csv"))
hr <- fread(file.path(dest, "pollutant_results.csv"))
stopifnot(nrow(hr) == 6L, nrow(c) == 63L,
  all(is.finite(hr$HR)), all(hr$HR_low <= hr$HR & hr$HR <= hr$HR_high),
  all(c$survival_without_transplant >= 0 & c$survival_without_transplant <= 1),
  all(c$survival_after_transplant >= 0 & c$survival_after_transplant <= 1),
  all(abs(c$benefit_probability - (c$survival_after_transplant - c$survival_without_transplant)) < 1e-10))
fwrite(rbindlist(checks), file.path(dest, "validation_checks.csv"))
fwrite(rbindlist(missing), file.path(dest, "complete_case_coverage.csv"))
cat("Real-data interval, mortality, exposure/LAS timing, model-event and paired-survival checks passed.\n")
