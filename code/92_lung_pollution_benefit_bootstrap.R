#!/usr/bin/env Rscript
# Resume uncertainty estimation/reporting using the completed model cache.
# The fitted cache is explicit, not silently reused by the main analysis script.
source("code/r_runtime.R"); ensure_user_library()
suppressPackageStartupMessages({library(data.table); library(survival); library(ggplot2)})
source("code/lung_survival_benefit_helpers.R")
setDTthreads(2)
dest <- Sys.getenv("LUNG_BENEFIT_OUTPUT", "output/lung_pollution_survival_benefit")
cache <- file.path(dest, "cache")
prepared <- readRDS(file.path(cache, "prepared.rds"))
fitted <- readRDS(file.path(cache, "fitted.rds"))
pre_d <- prepared$pre; post_d <- prepared$post
models <- fitted$models; refs <- fitted$reference; grids <- fitted$grids; point <- fitted$point
rm(prepared, fitted); gc(FALSE)
pollutants <- c("pm25", "no2", "o3")
horizon <- unique(point$horizon_days)
stopifnot(length(horizon) == 1L)
bootstrap_n <- as.integer(Sys.getenv("LUNG_BENEFIT_BOOTSTRAPS", "200"))
seed <- as.integer(Sys.getenv("LUNG_BENEFIT_SEED", "20261102"))
stopifnot(is.finite(bootstrap_n), bootstrap_n >= 0L)
log_msg <- function(...) message(format(Sys.time(), "%F %T"), " | ", paste0(..., collapse = ""))
# Reuse the main script's exact bootstrap/report statements, without executing
# preparation or point-model fitting. Stop before its methods/session reporting.
statements <- as.list(parse("code/91_lung_pollution_survival_benefit.R"))
is_call <- function(x, name) is.call(x) && identical(x[[1]], as.name(name))
start <- which(vapply(statements, function(x) is_call(x, "set.seed"), logical(1)))
end <- which(vapply(statements, function(x) is_call(x, "writeLines"), logical(1)))
stopifnot(length(start) == 1L, length(end) == 1L, end > start)
for (statement in statements[seq.int(start, end - 1L)]) eval(statement)
writeLines(c(paste("Requested person bootstrap resamples:", bootstrap_n), paste("Seed:", seed),
  "Completed against explicit fitted.rds cache; original point-model methods retained in methods.txt.",
  "Intervals condition on fixed empirical reference recipients. Inspect bootstrap_log.csv for failures."),
  file.path(dest, "bootstrap_run.txt"))
print(contrast)
