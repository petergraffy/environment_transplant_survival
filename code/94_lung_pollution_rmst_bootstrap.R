#!/usr/bin/env Rscript
# Joint person bootstrap of the two Cox models, with RMST integrated in EVERY
# refit. Parallel workers and resumable checkpoints keep the full-data run usable.
source("code/r_runtime.R"); ensure_user_library()
suppressPackageStartupMessages({library(data.table); library(survival); library(parallel)})
source("code/lung_survival_benefit_helpers.R")
setDTthreads(1)
dest <- Sys.getenv("LUNG_BENEFIT_OUTPUT", "output/lung_pollution_survival_benefit")
bootstrap_n <- as.integer(Sys.getenv("LUNG_RMST_BOOTSTRAPS", "200"))
seed <- as.integer(Sys.getenv("LUNG_BENEFIT_SEED", "20261102"))
workers <- as.integer(Sys.getenv("LUNG_RMST_WORKERS", as.character(min(12L, max(1L, detectCores() %/% 2L)))))
stopifnot(is.finite(bootstrap_n), bootstrap_n > 0L, is.finite(workers), workers > 0L,
  is.finite(seed), seed > 0, seed + bootstrap_n < .Machine$integer.max)
fitted_path <- normalizePath(file.path(dest, "cache", "fitted.rds"), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(dest, "rmst_1yr_exposure_curves.csv"))) stop("Run code/93_lung_pollution_rmst.R first")
# Checkpoints depend on source model, method code, package versions and seed.
# B and worker count are deliberately excluded: a run can be extended/resumed.
signature <- tempfile()
saveRDS(list(version = 1L, fitted_path = fitted_path, bytes = file.info(fitted_path)$size,
  modified = as.numeric(file.info(fitted_path)$mtime), seed = seed,
  code = tools::md5sum(c("code/lung_survival_benefit_helpers.R", "code/94_lung_pollution_rmst_bootstrap.R")),
  survival = as.character(packageVersion("survival")), R = R.version.string), signature)
key <- unname(tools::md5sum(signature)); unlink(signature)
checkpoint_dir <- file.path(dest, "cache", "rmst_bootstrap", key)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- normalizePath(checkpoint_dir, winslash = "/", mustWork = TRUE)
bootstrap_model_path <- file.path(checkpoint_dir, "bootstrap_models.rds")
if (!file.exists(bootstrap_model_path)) {
  compact <- readRDS(fitted_path)
  for (p in names(compact$models)) for (phase in names(compact$models[[p]])) {
    m <- compact$models[[p]][[phase]]
    compact$models[[p]][[phase]] <- list(data = m$data, fit = list(coefficients = coef(m$fit)))
  }
  saveRDS(compact, bootstrap_model_path, compress = FALSE)
  rm(compact, m); gc(FALSE)
}
job_file <- function(i) file.path(checkpoint_dir, sprintf("replicate_%04d.rds", i))
pending <- which(!file.exists(vapply(seq_len(bootstrap_n), job_file, character(1))))
message(format(Sys.time(), "%F %T"), " | RMST bootstrap: ", bootstrap_n,
  " replicates; ", length(pending), " pending; ", min(workers, length(pending)), " workers")
fwrite(data.table(seed, requested_resamples = bootstrap_n, workers,
  checkpoint_key = key, pending = length(pending)), file.path(dest, "rmst_bootstrap_run.csv"))

rmst_bootstrap_job <- function(iteration) {
  started <- proc.time()["elapsed"]
  set.seed(seed + iteration)
  counts <- tabulate(sample.int(length(bootstrap_ids), length(bootstrap_ids), replace = TRUE), nbins = length(bootstrap_ids))
  weights <- setNames(counts, as.character(bootstrap_ids))
  exposure <- list(); time <- list(); logs <- list()
  for (p in c("pm25", "no2", "o3")) {
    ans <- tryCatch({
      original_a <- fitted$models[[p]]$pretransplant
      original_b <- fitted$models[[p]]$posttransplant
      a <- lung_fit_cox(original_a$data, p, weights, coef(original_a$fit))
      b <- lung_fit_cox(original_b$data, p, weights, coef(original_b$fit))
      warnings <- c(a$warnings, b$warnings)
      if (any(grepl("converg|infinite", warnings, ignore.case = TRUE))) stop(paste(warnings, collapse = "; "))
      r <- fitted$reference[[p]]
      supported <- lung_supported_reference(r, a, b, 365, min_at_risk = 1L)
      if (nrow(supported) != nrow(r)) stop("Bootstrap lacks follow-up support for fixed reference population")
      h <- list(pre = lung_baseline(a), post = lung_baseline(b))
      ec <- lung_predict_rmst(a, b, r, p, fitted$grids[[p]], 365, h)
      tc <- lung_predict_rmst(a, b, r, p, range(fitted$grids[[p]]),
        sort(unique(c(seq(0, 360, by = 30), 365))), h)
      ec[, iteration := iteration]; tc[, iteration := iteration]
      list(exposure = ec, time = tc, warnings = paste(warnings, collapse = "; "))
    }, error = function(e) e)
    if (inherits(ans, "error")) {
      logs[[p]] <- data.table(iteration, pollutant = p, status = conditionMessage(ans), warnings = "")
    } else {
      logs[[p]] <- data.table(iteration, pollutant = p, status = "OK", warnings = ans$warnings)
      exposure[[p]] <- ans$exposure; time[[p]] <- ans$time
    }
    message("Replicate ", iteration, " / ", p, ": ", logs[[p]]$status)
    rm(ans); gc(FALSE)
  }
  result <- list(exposure = rbindlist(exposure), time = rbindlist(time), log = rbindlist(logs),
    elapsed_seconds = unname(proc.time()["elapsed"] - started))
  path <- file.path(checkpoint_dir, sprintf("replicate_%04d.rds", iteration))
  tmp <- paste0(path, ".tmp")
  saveRDS(result, tmp)
  if (!file.rename(tmp, path)) stop("Could not commit checkpoint ", path)
  iteration
}

if (length(pending)) {
  cl <- makePSOCKcluster(min(workers, length(pending)), outfile = file.path(dest, "cache", "rmst_bootstrap_workers.log"))
  tryCatch({
    clusterExport(cl, c("bootstrap_model_path", "seed", "checkpoint_dir", "rmst_bootstrap_job"), envir = globalenv())
    clusterEvalQ(cl, {
      source("code/r_runtime.R"); ensure_user_library()
      suppressPackageStartupMessages({library(data.table); library(survival)})
      source("code/lung_survival_benefit_helpers.R")
      setDTthreads(1)
      fitted <- readRDS(bootstrap_model_path)
      bootstrap_ids <- sort(unique(unlist(lapply(fitted$models, function(m)
        c(m$pretransplant$data$PERS_ID, m$posttransplant$data$PERS_ID)))))
      NULL
    })
    invisible(parLapplyLB(cl, pending, rmst_bootstrap_job))
  }, finally = stopCluster(cl))
}

replicates <- lapply(seq_len(bootstrap_n), function(i) readRDS(job_file(i)))
logs <- rbindlist(lapply(replicates, `[[`, "log"))
boot_exposure <- rbindlist(lapply(replicates, `[[`, "exposure"))
boot_time <- rbindlist(lapply(replicates, `[[`, "time"))
fwrite(logs, file.path(dest, "rmst_bootstrap_log.csv"))
fwrite(boot_exposure, file.path(dest, "rmst_bootstrap_exposure_curves.csv"))
fwrite(boot_time, file.path(dest, "rmst_bootstrap_time_curves.csv"))
fwrite(data.table(iteration = seq_len(bootstrap_n), elapsed_seconds = vapply(replicates, `[[`, numeric(1), "elapsed_seconds")),
  file.path(dest, "rmst_bootstrap_timing.csv"))

attach_intervals <- function(point, boot) {
  point <- copy(point)
  old <- grep("(_low|_high|_se)$|^successful_bootstraps$|^confidence_status$", names(point), value = TRUE)
  if (length(old)) point[, (old) := NULL]
  keys <- c("pollutant", "concentration", "horizon_days")
  metrics <- c("rmst_without_transplant_days", "rmst_after_transplant_days", "rmst_benefit_days")
  if (nrow(boot)) {
    bands <- boot[, c(list(successful_bootstraps = uniqueN(iteration)),
      unlist(lapply(.SD, function(x) list(low = unname(quantile(x, .025)), high = unname(quantile(x, .975)), se = sd(x))), recursive = FALSE)),
      by = keys, .SDcols = metrics]
    setnames(bands, setdiff(names(bands), c(keys, "successful_bootstraps")),
      as.vector(t(outer(metrics, c("low", "high", "se"), paste, sep = "_"))))
    point <- merge(point, bands, by = keys, all.x = TRUE)
    point[is.na(successful_bootstraps), successful_bootstraps := 0L]
  } else point[, successful_bootstraps := 0L]
  point[, confidence_status := fifelse(successful_bootstraps >= max(100L, ceiling(.8 * bootstrap_n)),
    "estimated", "insufficient_resamples")]
  for (nm in grep("_(low|high|se)$", names(point), value = TRUE)) point[confidence_status != "estimated", (nm) := NA_real_]
  point[]
}
curves <- attach_intervals(fread(file.path(dest, "rmst_1yr_exposure_curves.csv")), boot_exposure)
times <- attach_intervals(fread(file.path(dest, "rmst_1yr_time_curves.csv")), boot_time)
fwrite(curves, file.path(dest, "rmst_1yr_exposure_curves.csv"))
fwrite(times, file.path(dest, "rmst_1yr_time_curves.csv"))
contrast <- curves[order(concentration), .(low_exposure = first(concentration), high_exposure = last(concentration),
  rmst_benefit_low_days = first(rmst_benefit_days), rmst_benefit_high_days = last(rmst_benefit_days),
  benefit_difference_days = last(rmst_benefit_days) - first(rmst_benefit_days), recipients = first(recipients),
  confidence_status = first(confidence_status)), by = pollutant]
if (nrow(boot_exposure)) {
  delta <- boot_exposure[order(concentration), .(difference = last(rmst_benefit_days) - first(rmst_benefit_days)), by = .(pollutant, iteration)]
  dc <- delta[, .(benefit_difference_days_low = unname(quantile(difference, .025)),
    benefit_difference_days_high = unname(quantile(difference, .975)), successful_bootstraps = .N), by = pollutant]
  contrast <- merge(contrast, dc, by = "pollutant", all.x = TRUE)
  contrast[confidence_status != "estimated", c("benefit_difference_days_low", "benefit_difference_days_high") := .(NA_real_, NA_real_)]
}
fwrite(contrast, file.path(dest, "rmst_1yr_exposure_contrasts.csv"))
# Rerender from aggregate tables, without rereading/refitting the clinical cache.
suppressPackageStartupMessages(library(ggplot2))
statements <- as.list(parse("code/93_lung_pollution_rmst.R"))
start <- which(vapply(statements, function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
  identical(x[[2]], as.name("labels")), logical(1)))
for (statement in statements[seq.int(start, length(statements))]) eval(statement)
writeLines(c("Pointwise 95% percentile intervals from joint person bootstrap with both Cox models refitted in every replicate.",
  "The same person multiplicity weights are used across phases and pollutants within each replicate; recipient reference populations remain fixed.",
  "RMST is reintegrated from the refitted baseline hazards and coefficients, including their joint uncertainty and cross-phase dependence.",
  "Covariates and exposure scenarios are frozen as in the original RMST analysis. Warm starts change computation, not the fitted model specification.",
  "Checkpoints permit exact resumption/extension. Replicate seeds equal base seed plus replicate number, independent of worker scheduling.",
  paste("Requested replicates:", bootstrap_n, "; base seed:", seed, "; checkpoint key:", key),
  "At least 100 and 80% of requested replicates must succeed per pollutant before confidence intervals are displayed.",
  "Bands are pointwise, not simultaneous, and conditional on the fixed empirical recipient reference population; they do not remove observational selection or confounding."),
  file.path(dest, "rmst_bootstrap_methods.txt"))
print(logs[, .N, by = .(pollutant, status)])
