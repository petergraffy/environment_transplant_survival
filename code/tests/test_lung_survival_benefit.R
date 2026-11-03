# Synthetic tests: no restricted source data required.
source("code/r_runtime.R"); ensure_user_library()
suppressPackageStartupMessages({library(data.table); library(survival)})
source("code/lung_survival_benefit_helpers.R")

cohort <- data.table(PERS_ID = 1:4,
  index_date = as.Date(c("2020-01-15", "2020-01-15", "2020-02-01", "2020-01-31")),
  end_date = as.Date(c("2020-03-01", "2020-01-15", "2020-02-10", "2020-02-02")),
  death_event = c(1L, 1L, 0L, 1L))
d <- lung_month_intervals(cohort)
stopifnot(all(d$tstop > d$tstart), sum(d$event) == 3L,
  d[PERS_ID == 2L, tstop] == .5,
  all(d[, .(contiguous = all(tstart[-1] == head(tstop, -1))), by = PERS_ID]$contiguous),
  d[PERS_ID == 1L, sum(tstop - tstart)] == 46,
  d[PERS_ID == 1L & event == 1L, start_date] == as.Date("2020-02-01"))
d[, candidate_zip := "00123"]
monthly <- data.table(zip = "00123", year = c(2019L, 2020L), month = c(12L, 1L),
  pm25 = c(4, 8), no2 = c(10, 20), o3 = c(30, 40), no2_annual_fallback = FALSE,
  exposure_source_end = as.Date(c("2019-12-31", "2020-01-31")))
e <- lung_attach_exposure(d, monthly)
stopifnot(all(e[format(start_date, "%m") == "01", pm25] == 4),
  all(e[format(start_date, "%m") == "02", pm25] == 8))
scores <- data.table(PERS_ID = c(1L, 1L, 1L), score_date = as.Date(c("2020-01-01", "2020-02-01", "2020-03-01")),
  las = c(30, 40, 90))
s <- lung_attach_las(d, scores, 35)
stopifnot(s[PERS_ID == 1L, identical(las_per10, c(3, 4))],
  all(s[PERS_ID != 1L, las_unavailable] == 1L),
  all(s[PERS_ID != 1L, las_per10] == 3.5))

# Compare our conditional-survival calculation with survival::survfit, including
# center strata, nonzero waiting time and distinct listing/transplant ages.
set.seed(91)
n <- 800L
x <- data.table(PERS_ID = seq_len(n), tstart = 0, tstop = pmin(rexp(n, 1 / 700), 1600),
  age = runif(n, 20, 70), sex = sample(c("F", "M"), n, TRUE),
  race = sample(c("A", "B"), n, TRUE), zcta_svi_proxy = runif(n),
  listing_center = sample(c("center_A", "center_B"), n, TRUE), las_per10 = runif(n, 3, 7),
  las_unavailable = 0L, pm25 = runif(n, 4, 16), no2_annual_fallback = FALSE)
x[, event := as.integer(tstop < 1600)]
a <- lung_fit_cox(x, "pm25")
z <- lung_fit_cox(copy(x)[, tstop := tstop * .85], "pm25")
r <- x[1:3]
r[, `:=`(listing_age = age, transplant_age = age + .5, waiting_days = c(50, 100, 200))]
q <- lung_predict_benefit(a, z, r, "pm25", c(6, 12))
stopifnot(nrow(q) == 2L, all(q$survival_without_transplant >= 0 & q$survival_without_transplant <= 1),
  all(abs(q$benefit_probability - (q$survival_after_transplant - q$survival_without_transplant)) < 1e-12))
aa <- as.data.frame(r[1]); aa$exposure_scaled <- 6 / 5
sf <- survfit(a$fit, newdata = aa)
sv <- summary(sf, times = c(50, 415), extend = FALSE)$surv
hh <- lung_baseline(a)
lp <- predict(a$fit, newdata = aa, type = "lp", reference = "zero")
ours <- exp(-(lung_hazard_at(hh, aa$listing_center, 415) - lung_hazard_at(hh, aa$listing_center, 50)) * exp(lp))
stopifnot(length(sv) == 2L, abs(ours - sv[2] / sv[1]) < 1e-10)
# Support requires follow-up to waiting time PLUS horizon, not just horizon.
unsupported <- copy(r); unsupported[, waiting_days := 10000]
stopifnot(nrow(lung_supported_reference(unsupported, a, z)) == 0L)
# Bootstrap refits must work with multiplicity weights and the same reference.
w <- setNames(tabulate(sample.int(n, n, TRUE), nbins = n), as.character(seq_len(n)))
ab <- lung_fit_cox(x, "pm25", w)
stopifnot(is.finite(coef(ab$fit)["exposure_scaled"]))
cat("Lung benefit tests passed: interval boundaries, perioperative death, lagged exposure, LAS timing, conditional survival, horizon support and person bootstrap.\n")

# Exercise all six fits plus the standalone bootstrap/report entry point using
# synthetic data, including successful-resample accounting and CI suppression
# for a smoke run too small to support reported confidence intervals.
synthetic_dest <- file.path("output/lung_pollution_survival_benefit/cache", "synthetic_validation")
dir.create(file.path(synthetic_dest, "cache"), recursive = TRUE, showWarnings = FALSE)
x[, `:=`(no2 = pm25 * 1.5, o3 = pm25 + 25)]
synthetic_models <- list(); synthetic_refs <- list(); synthetic_grids <- list(); synthetic_point <- list()
for (p in c("pm25", "no2", "o3")) {
  aa <- lung_fit_cox(x, p)
  zz <- lung_fit_cox(copy(x)[, tstop := tstop * .85], p)
  synthetic_models[[p]] <- list(pretransplant = aa, posttransplant = zz)
  synthetic_refs[[p]] <- r
  synthetic_grids[[p]] <- range(x[[p]])
  synthetic_point[[p]] <- lung_predict_benefit(aa, zz, r, p, synthetic_grids[[p]])
}
saveRDS(list(pre = x, post = x, reference = r), file.path(synthetic_dest, "cache", "prepared.rds"))
saveRDS(list(models = synthetic_models, reference = synthetic_refs, grids = synthetic_grids,
  point = rbindlist(synthetic_point)), file.path(synthetic_dest, "cache", "fitted.rds"))
old_output <- Sys.getenv("LUNG_BENEFIT_OUTPUT", unset = NA_character_)
old_boot <- Sys.getenv("LUNG_BENEFIT_BOOTSTRAPS", unset = NA_character_)
Sys.setenv(LUNG_BENEFIT_OUTPUT = synthetic_dest, LUNG_BENEFIT_BOOTSTRAPS = "5")
runner <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
status <- system2(runner, "code/92_lung_pollution_benefit_bootstrap.R",
  stdout = file.path(synthetic_dest, "run.log"), stderr = file.path(synthetic_dest, "run_errors.log"))
if (is.na(old_output)) Sys.unsetenv("LUNG_BENEFIT_OUTPUT") else Sys.setenv(LUNG_BENEFIT_OUTPUT = old_output)
if (is.na(old_boot)) Sys.unsetenv("LUNG_BENEFIT_BOOTSTRAPS") else Sys.setenv(LUNG_BENEFIT_BOOTSTRAPS = old_boot)
if (status != 0L) stop(paste(readLines(file.path(synthetic_dest, "run_errors.log")), collapse = "\n"))
out <- fread(file.path(synthetic_dest, "survival_benefit_curves.csv"))
bl <- fread(file.path(synthetic_dest, "bootstrap_log.csv"))
stopifnot(nrow(bl) == 15L, all(bl$status == "OK"), all(out$successful_bootstraps == 5L),
  all(out$confidence_status == "insufficient_resamples"), all(is.na(out$benefit_probability_low)))
cat("Standalone six-model bootstrap/report smoke test passed.\n")

# If one pollutant has no successful resamples, retain its point curve and
# explicitly mark missing uncertainty instead of silently dropping that panel.
sfitted <- readRDS(file.path(synthetic_dest, "cache", "fitted.rds"))
sfitted$models$o3$pretransplant$data[, event := 0L]
saveRDS(sfitted, file.path(synthetic_dest, "cache", "fitted.rds"))
Sys.setenv(LUNG_BENEFIT_OUTPUT = synthetic_dest, LUNG_BENEFIT_BOOTSTRAPS = "2")
status <- system2(runner, "code/92_lung_pollution_benefit_bootstrap.R",
  stdout = file.path(synthetic_dest, "failure_test.log"), stderr = file.path(synthetic_dest, "failure_test_errors.log"))
if (is.na(old_output)) Sys.unsetenv("LUNG_BENEFIT_OUTPUT") else Sys.setenv(LUNG_BENEFIT_OUTPUT = old_output)
if (is.na(old_boot)) Sys.unsetenv("LUNG_BENEFIT_BOOTSTRAPS") else Sys.setenv(LUNG_BENEFIT_BOOTSTRAPS = old_boot)
stopifnot(status == 0L)
out <- fread(file.path(synthetic_dest, "survival_benefit_curves.csv"))
stopifnot(nrow(out) == 6L, all(out[pollutant == "o3", successful_bootstraps] == 0L),
  all(out[pollutant == "o3", confidence_status] == "not_run_or_all_failed"))
cat("Failed-bootstrap point-curve preservation test passed.\n")
