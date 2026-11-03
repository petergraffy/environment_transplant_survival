# Validate the completed clinical-data RMST confidence intervals independently
# against saved replicate estimates. This does not refit clinical models.
source("code/r_runtime.R"); ensure_user_library()
suppressPackageStartupMessages(library(data.table))
dest <- Sys.getenv("LUNG_BENEFIT_OUTPUT", "output/lung_pollution_survival_benefit")
curves <- fread(file.path(dest, "rmst_1yr_exposure_curves.csv"))
times <- fread(file.path(dest, "rmst_1yr_time_curves.csv"))
boot <- fread(file.path(dest, "rmst_bootstrap_exposure_curves.csv"))
boot_time <- fread(file.path(dest, "rmst_bootstrap_time_curves.csv"))
logs <- fread(file.path(dest, "rmst_bootstrap_log.csv"))
run <- fread(file.path(dest, "rmst_bootstrap_run.csv"))
contrasts <- fread(file.path(dest, "rmst_1yr_exposure_contrasts.csv"))
stopifnot(nrow(curves) == 63L, nrow(times) == 84L, all(curves$confidence_status == "estimated"),
  all(times$confidence_status == "estimated"), all(curves$successful_bootstraps >= max(100, ceiling(.8 * run$requested_resamples))))
metrics <- c("rmst_without_transplant_days", "rmst_after_transplant_days", "rmst_benefit_days")
keys <- c("pollutant", "concentration", "horizon_days")
for (pair in list(list(point = curves, boot = boot), list(point = times, boot = boot_time))) {
  x <- pair$point; b <- pair$boot
  stopifnot(all(abs(b$rmst_benefit_days - (b$rmst_after_transplant_days - b$rmst_without_transplant_days)) < 1e-8))
  for (metric in metrics) {
    expected <- b[, .(low = unname(quantile(get(metric), .025)), high = unname(quantile(get(metric), .975)),
      se = sd(get(metric)), n = uniqueN(iteration)), by = keys]
    actual <- x[, c(keys, paste0(metric, c("_low", "_high", "_se")), "successful_bootstraps"), with = FALSE]
    setnames(actual, paste0(metric, c("_low", "_high", "_se")), c("actual_low", "actual_high", "actual_se"))
    v <- merge(expected, actual, by = keys)
    stopifnot(nrow(v) == nrow(x), all(is.finite(v$actual_low)), all(is.finite(v$actual_high)),
      all(v$actual_low <= v$actual_high), max(abs(v$low - v$actual_low)) < 1e-8,
      max(abs(v$high - v$actual_high)) < 1e-8, max(abs(v$se - v$actual_se)) < 1e-8,
      all(v$n == v$successful_bootstraps))
  }
  for (metric in metrics[1:2]) stopifnot(all(x[[paste0(metric, "_low")]] >= -1e-8),
    all(x[[paste0(metric, "_high")]] <= x$horizon_days + 1e-8))
}
zero <- times[horizon_days == 0]
stopifnot(all(abs(unlist(zero[, c(metrics, paste0(metrics, "_low"), paste0(metrics, "_high")), with = FALSE])) < 1e-10))
delta <- boot[order(concentration), .(difference = last(rmst_benefit_days) - first(rmst_benefit_days)), by = .(pollutant, iteration)]
dc <- delta[, .(expected_low = unname(quantile(difference, .025)), expected_high = unname(quantile(difference, .975))), by = pollutant]
v <- merge(contrasts, dc, by = "pollutant")
stopifnot(nrow(v) == 3L, max(abs(v$benefit_difference_days_low - v$expected_low)) < 1e-8,
  max(abs(v$benefit_difference_days_high - v$expected_high)) < 1e-8)
counts <- logs[status == "OK", .(successful = .N), by = pollutant]
cc <- merge(unique(curves[, .(pollutant, successful_bootstraps)]), counts, by = "pollutant")
stopifnot(all(cc$successful_bootstraps == cc$successful))
figures <- as.vector(outer(c("benefit", "components", "time"), c("png", "pdf"),
  function(x, y) paste0("rmst_1yr_", x, ".", y)))
stopifnot(all(file.exists(file.path(dest, figures))))
fwrite(data.table(check = c("All exposure/time curves have adequate bootstrap success", "Reported percentile limits and SEs match saved replicates",
  "Replicate benefits equal paired RMST differences", "RMST intervals obey physical bounds and zero-time identity",
  "Exposure contrast limits use within-replicate paired differences", "Success counts match replicate logs", "All six figure files exist"), passed = TRUE),
  file.path(dest, "rmst_confidence_validation.csv"))
print(cc)
cat("Clinical RMST confidence interval validation passed.\n")
