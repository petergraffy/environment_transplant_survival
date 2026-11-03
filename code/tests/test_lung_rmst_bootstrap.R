source("code/r_runtime.R"); ensure_user_library()
suppressPackageStartupMessages({library(data.table); library(survival)})
source("code/lung_survival_benefit_helpers.R")
set.seed(94)
n <- 240L
d <- data.table(PERS_ID = seq_len(n), tstart = 0, tstop = runif(n, 10, 1000), event = sample(0:1, n, TRUE),
  pm25 = runif(n, 4, 12), age = runif(n, 20, 70), sex = sample(c("F", "M"), n, TRUE),
  race = sample(c("A", "B"), n, TRUE), zcta_svi_proxy = runif(n), listing_center = sample(c("a", "b"), n, TRUE),
  las_per10 = runif(n, 3, 6), las_unavailable = 0L)
d[, `:=`(no2 = pm25 * 1.5, o3 = pm25 + 25)]
models <- list(); reference <- list(); grids <- list(); points <- list()
for (p in c("pm25", "no2", "o3")) {
  m <- lung_fit_cox(d, p)
  warm <- lung_fit_cox(d, p, setNames(rep(1, n), as.character(d$PERS_ID)), coef(m$fit))
  stopifnot(max(abs(coef(m$fit) - coef(warm$fit))) < 1e-6)
  models[[p]] <- list(pretransplant = m, posttransplant = m)
  r <- copy(d[1:5]); r[, `:=`(listing_age = age, transplant_age = age, waiting_days = 0)]
  reference[[p]] <- r; grids[[p]] <- range(d[[p]])
  points[[p]] <- lung_predict_benefit(m, m, r, p, grids[[p]])
}
test_dest <- "output/lung_pollution_survival_benefit/cache/rmst_bootstrap_validation"
dir.create(file.path(test_dest, "cache"), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(models = models, reference = reference, grids = grids, point = rbindlist(points)),
  file.path(test_dest, "cache", "fitted.rds"))
old <- Sys.getenv(c("LUNG_BENEFIT_OUTPUT", "LUNG_RMST_BOOTSTRAPS", "LUNG_RMST_WORKERS"), unset = NA_character_)
Sys.setenv(LUNG_BENEFIT_OUTPUT = test_dest, LUNG_RMST_BOOTSTRAPS = "2", LUNG_RMST_WORKERS = "2")
runner <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
run_script <- function(path, label) {
  status <- system2(runner, path, stdout = file.path(test_dest, paste0(label, ".log")),
    stderr = file.path(test_dest, paste0(label, "_errors.log")))
  if (status != 0L) stop(paste(readLines(file.path(test_dest, paste0(label, "_errors.log"))), collapse = "\n"))
}
run_script("code/93_lung_pollution_rmst.R", "point")
run_script("code/94_lung_pollution_rmst_bootstrap.R", "bootstrap")
boot <- fread(file.path(test_dest, "rmst_bootstrap_exposure_curves.csv"))
log <- fread(file.path(test_dest, "rmst_bootstrap_log.csv"))
stopifnot(nrow(log) == 6L, all(log$status == "OK"), all(abs(boot$rmst_benefit_days) < 1e-9))
key <- fread(file.path(test_dest, "rmst_bootstrap_run.csv"))$checkpoint_key
files <- list.files(file.path(test_dest, "cache", "rmst_bootstrap", key), "[.]rds$", full.names = TRUE)
mtime <- file.info(files)$mtime
run_script("code/94_lung_pollution_rmst_bootstrap.R", "resume")
stopifnot(fread(file.path(test_dest, "rmst_bootstrap_run.csv"))$pending == 0L,
  identical(mtime, file.info(files)$mtime))
for (name in names(old)) if (is.na(old[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, setNames(list(old[[name]]), name))

# Independently check percentile extraction, SEs and insufficient-resample
# protection using a known set of 200 estimates, without clinical refits.
statements <- as.list(parse("code/94_lung_pollution_rmst_bootstrap.R"))
fun <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
  identical(x[[2]], as.name("attach_intervals")), statements)
stopifnot(length(fun) == 1L)
eval(fun[[1]])
bootstrap_n <- 200L
point <- data.table(pollutant = "pm25", concentration = 5, horizon_days = 365,
  rmst_without_transplant_days = 320, rmst_after_transplant_days = 338, rmst_benefit_days = 18)
b <- point[rep(1, 200)]
b[, `:=`(iteration = seq_len(.N), rmst_without_transplant_days = 320 + seq_len(.N) / 100,
  rmst_after_transplant_days = 338 + seq_len(.N) / 50, rmst_benefit_days = 18 + seq_len(.N) / 100)]
ci <- attach_intervals(point, b)
stopifnot(ci$confidence_status == "estimated", ci$successful_bootstraps == 200L,
  abs(ci$rmst_benefit_days_low - quantile(b$rmst_benefit_days, .025)) < 1e-10,
  abs(ci$rmst_benefit_days_high - quantile(b$rmst_benefit_days, .975)) < 1e-10,
  abs(ci$rmst_benefit_days_se - sd(b$rmst_benefit_days)) < 1e-10,
  all(is.na(attach_intervals(point, b[1:5])$rmst_benefit_days_low)))
suppressPackageStartupMessages(library(ggplot2))
dest <- test_dest
curves <- fread(file.path(dest, "rmst_1yr_exposure_curves.csv"))
times <- fread(file.path(dest, "rmst_1yr_time_curves.csv"))
contrast <- fread(file.path(dest, "rmst_1yr_exposure_contrasts.csv"))
for (table in list(curves, times)) {
  table[, confidence_status := "estimated"]
  for (metric in c("rmst_without_transplant_days", "rmst_after_transplant_days", "rmst_benefit_days")) {
    table[, (paste0(metric, "_low")) := get(metric) - pmin(.1, horizon_days / 365)]
    table[, (paste0(metric, "_high")) := get(metric) + pmin(.1, horizon_days / 365)]
  }
}
plot_statements <- as.list(parse("code/93_lung_pollution_rmst.R"))
start <- which(vapply(plot_statements, function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
  identical(x[[2]], as.name("labels")), logical(1)))
for (statement in plot_statements[seq.int(start, length(plot_statements))]) eval(statement)
stopifnot(has_ci, all(c("ymin", "ymax") %in% names(ggplot_build(benefit)$data[[1]])),
  all(c("ymin", "ymax") %in% names(ggplot_build(component_plot)$data[[1]])),
  all(c("ymin", "ymax") %in% names(ggplot_build(time_plot)$data[[1]])))
cat("RMST bootstrap tests passed: warm starts, paired person resampling, parallel execution, checkpoint resumption and 95% percentile intervals.\n")
