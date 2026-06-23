#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(readr)
})

pm25_dir <- file.path("data", "processed", "lghap_pm25_zcta_daily")
aod_dir <- file.path("data", "processed", "lghap_aod_zcta_daily")
out_dir <- file.path("data", "processed", "lghap_pm25_zcta_daily_derived")
fig_dir <- file.path("output", "figures", "lghap_pm25_derived_2022")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

calibration_year <- as.integer(Sys.getenv("LGHAP_DERIVED_PM25_CALIBRATION_YEAR", "2021"))
target_year <- as.integer(Sys.getenv("LGHAP_DERIVED_PM25_TARGET_YEAR", "2022"))
max_train_rows <- as.integer(Sys.getenv("LGHAP_DERIVED_PM25_MAX_TRAIN_ROWS", "2000000"))
seed <- as.integer(Sys.getenv("LGHAP_DERIVED_PM25_SEED", "20250605"))

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

read_dt <- function(path, select = NULL) {
  if (grepl("\\.gz$", path) && !requireNamespace("R.utils", quietly = TRUE)) {
    dt <- as.data.table(read_csv(path, show_col_types = FALSE, progress = FALSE))
    if (!is.null(select)) dt <- dt[, ..select]
    return(dt)
  }
  fread(path, select = select)
}

read_monthly <- function(dir, prefix, year, months = 1:12) {
  files <- file.path(dir, sprintf("%s_%04d_%02d.csv.gz", prefix, year, months))
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Missing required monthly files:\n", paste(missing, collapse = "\n"), call. = FALSE)
  rbindlist(lapply(files, read_dt), fill = TRUE)
}

log_msg("Reading historical PM2.5 for monthly ZCTA climatology")
hist_files <- list.files(
  pm25_dir,
  pattern = "^lghap_pm25_zcta_daily_20(0[5-9]|1[0-9]|20|21)_[0-9]{2}\\.csv\\.gz$",
  full.names = TRUE
)
if (!length(hist_files)) stop("No historical LGHAP PM2.5 daily ZCTA files found in ", pm25_dir, call. = FALSE)
clim_parts <- lapply(seq_along(hist_files), function(i) {
  path <- hist_files[i]
  if (i %% 12L == 1L) log_msg("  climatology file ", i, "/", length(hist_files), " ", basename(path))
  dt <- read_dt(path, select = c("zip", "month", "pm25_ug_m3"))
  dt[, zip := sprintf("%05s", zip)]
  dt[, .(sum_pm25 = sum(pm25_ug_m3, na.rm = TRUE), n_pm25 = sum(is.finite(pm25_ug_m3))), by = .(zip, month)]
})
clim <- rbindlist(clim_parts, fill = TRUE)[
  ,
  .(pm25_clim_month = sum(sum_pm25, na.rm = TRUE) / sum(n_pm25, na.rm = TRUE)),
  by = .(zip, month)
]
rm(clim_parts)
gc()
fwrite(clim, file.path(out_dir, "lghap_pm25_zcta_monthly_climatology_2005_2021.csv"))

log_msg("Reading calibration year paired AOD and PM2.5: ", calibration_year)
pm_cal <- read_monthly(pm25_dir, "lghap_pm25_zcta_daily", calibration_year)
pm_cal <- pm_cal[, .(zip = sprintf("%05s", zip), date = as.Date(date), pm25_ug_m3)]
aod_cal <- read_monthly(aod_dir, "lghap_aod_zcta_daily", calibration_year)
aod_cal <- aod_cal[, .(zip = sprintf("%05s", zip), date = as.Date(date), aod)]
train <- merge(pm_cal, aod_cal, by = c("zip", "date"), all = FALSE)
train[, month := as.integer(format(date, "%m"))]
train <- merge(train, clim, by = c("zip", "month"), all.x = TRUE, sort = FALSE)
train <- train[is.finite(pm25_ug_m3) & is.finite(aod) & is.finite(pm25_clim_month)]
train[, `:=`(aod2 = aod^2, log1p_aod = log1p(pmax(aod, 0)))]

set.seed(seed)
train_fit <- train
if (max_train_rows > 0L && nrow(train_fit) > max_train_rows) {
  log_msg("Sampling ", max_train_rows, " rows for model fit from ", nrow(train_fit), " paired rows")
  train_fit <- train_fit[sample.int(nrow(train_fit), max_train_rows)]
}
train_fit[, fold := sample.int(5L, .N, replace = TRUE)]
train[, fold := sample.int(5L, .N, replace = TRUE)]

model_formula <- pm25_ug_m3 ~ pm25_clim_month + aod + aod2 + log1p_aod + factor(month) + pm25_clim_month:aod
log_msg("Fitting AOD-to-PM2.5 calibration model")
fit <- lm(model_formula, data = train_fit)
saveRDS(fit, file.path(out_dir, sprintf("lghap_pm25_from_aod_calibration_model_%04d.rds", calibration_year)))

coef_dt <- data.table(term = names(coef(fit)), estimate = as.numeric(coef(fit)))
fwrite(coef_dt, file.path(out_dir, sprintf("lghap_pm25_from_aod_calibration_coefficients_%04d.csv", calibration_year)))

log_msg("Running 5-fold calibration QA")
fold_metrics <- rbindlist(lapply(1:5, function(k) {
  fold_train <- train_fit[fold != k]
  fold_test <- train[fold == k]
  fold_fit <- lm(model_formula, data = fold_train)
  pred <- pmax(0, as.numeric(predict(fold_fit, newdata = fold_test)))
  err <- pred - fold_test$pm25_ug_m3
  data.table(
    fold = k,
    n = length(err),
    rmse = sqrt(mean(err^2, na.rm = TRUE)),
    mae = mean(abs(err), na.rm = TRUE),
    bias = mean(err, na.rm = TRUE),
    r2 = 1 - sum(err^2, na.rm = TRUE) / sum((fold_test$pm25_ug_m3 - mean(fold_test$pm25_ug_m3, na.rm = TRUE))^2, na.rm = TRUE)
  )
}))
overall_pred <- pmax(0, as.numeric(predict(fit, newdata = train)))
overall_err <- overall_pred - train$pm25_ug_m3
overall_metrics <- data.table(
  fold = 0L,
  n = length(overall_err),
  rmse = sqrt(mean(overall_err^2, na.rm = TRUE)),
  mae = mean(abs(overall_err), na.rm = TRUE),
  bias = mean(overall_err, na.rm = TRUE),
  r2 = 1 - sum(overall_err^2, na.rm = TRUE) / sum((train$pm25_ug_m3 - mean(train$pm25_ug_m3, na.rm = TRUE))^2, na.rm = TRUE)
)
metrics <- rbind(overall_metrics, fold_metrics)
fwrite(metrics, file.path(out_dir, sprintf("lghap_pm25_from_aod_calibration_qa_%04d.csv", calibration_year)))

qa_plot_dt <- train[sample.int(nrow(train), min(nrow(train), 250000))]
qa_plot_dt[, predicted_pm25_ug_m3 := pmax(0, as.numeric(predict(fit, newdata = qa_plot_dt)))]
ggsave(
  file.path(fig_dir, sprintf("lghap_pm25_from_aod_calibration_scatter_%04d.png", calibration_year)),
  ggplot(qa_plot_dt, aes(pm25_ug_m3, predicted_pm25_ug_m3)) +
    geom_bin2d(bins = 120) +
    geom_abline(slope = 1, intercept = 0, color = "white", linewidth = 0.6) +
    scale_fill_viridis_c(option = "magma", trans = "log10") +
    labs(
      title = sprintf("LGHAP AOD-Derived PM2.5 Calibration QA, %04d", calibration_year),
      x = "True LGHAP PM2.5 (ug/m3)",
      y = "Model-derived PM2.5 (ug/m3)",
      fill = "Rows"
    ) +
    theme_minimal(base_size = 12),
  width = 8,
  height = 6,
  dpi = 220,
  bg = "white"
)

log_msg("Predicting target year PM2.5 from AOD: ", target_year)
manifest <- rbindlist(lapply(1:12, function(month) {
  aod_path <- file.path(aod_dir, sprintf("lghap_aod_zcta_daily_%04d_%02d.csv.gz", target_year, month))
  if (!file.exists(aod_path)) stop("Missing target AOD monthly file: ", aod_path, call. = FALSE)
  pred_path <- file.path(out_dir, sprintf("lghap_pm25_zcta_daily_derived_from_aod_%04d_%02d.csv.gz", target_year, month))

  log_msg("  predicting ", target_year, "-", sprintf("%02d", month))
  dt <- read_dt(aod_path)
  dt <- dt[, .(zip = sprintf("%05s", zip), date = as.Date(date), year, month, aod, aod_value_source = value_source)]
  dt <- merge(dt, clim, by = c("zip", "month"), all.x = TRUE, sort = FALSE)
  dt[, `:=`(aod2 = aod^2, log1p_aod = log1p(pmax(aod, 0)))]
  dt[, pm25_ug_m3 := pmax(0, as.numeric(predict(fit, newdata = dt)))]
  dt[, `:=`(
    pollutant = "pm25_lghap_daily_derived_from_aod",
    temporal_resolution = "daily",
    derivation = sprintf("AOD-to-PM2.5 model calibrated on %04d LGHAP AOD and PM2.5 ZCTA aggregates", calibration_year),
    calibration_year = calibration_year,
    value_source = fifelse(is.na(pm25_ug_m3), "missing", "model_derived_from_aod")
  )]
  out <- dt[, .(
    zip, date, pollutant, year, month, temporal_resolution,
    pm25_ug_m3, value_source, aod, aod_value_source,
    pm25_clim_month, derivation, calibration_year
  )]
  setorder(out, zip, date)
  fwrite(out, pred_path)

  data.table(
    year = target_year,
    month = month,
    rows = nrow(out),
    zips = uniqueN(out$zip),
    dates = uniqueN(out$date),
    missing_values = sum(is.na(out$pm25_ug_m3)),
    mean = mean(out$pm25_ug_m3, na.rm = TRUE),
    p01 = quantile(out$pm25_ug_m3, 0.01, na.rm = TRUE),
    p50 = quantile(out$pm25_ug_m3, 0.50, na.rm = TRUE),
    p99 = quantile(out$pm25_ug_m3, 0.99, na.rm = TRUE),
    output = normalizePath(pred_path, winslash = "/", mustWork = TRUE)
  )
}), fill = TRUE)

fwrite(manifest, file.path(out_dir, sprintf("lghap_pm25_zcta_daily_derived_from_aod_manifest_%04d.csv", target_year)))
fwrite(manifest, file.path(out_dir, sprintf("lghap_pm25_zcta_daily_derived_from_aod_qc_%04d_all_months.csv", target_year)))
log_msg("Done. Derived PM2.5 outputs written to ", normalizePath(out_dir, winslash = "/"))
