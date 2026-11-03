source("code/r_runtime.R")
ensure_user_library()
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(survival)
  library(ggplot2); library(tibble); library(broom)
})
revision_dir <- "output/revision_2026_checklist"
dir.create(revision_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(revision_dir, "cache"), showWarnings = FALSE)
log_revision <- function(...) message(format(Sys.time(), "%F %T"), " | ", paste0(..., collapse = ""))

# Reuse the existing preparation code without executing its model/output section.
load_script_prefix <- function(path, stop_assignment, env = new.env(parent = globalenv())) {
  statements <- as.list(parse(path))
  boundary <- which(vapply(statements, function(x) is.call(x) &&
    identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name(stop_assignment)), logical(1)))
  stopifnot(length(boundary) >= 1L)
  boundary <- boundary[1L]
  stopifnot(boundary > 1L)
  for (statement in statements[seq_len(boundary - 1L)]) eval(statement, env)
  env
}

baseline_data <- function() {
  cache <- file.path(revision_dir, "cache", "baseline.rds")
  if (file.exists(cache)) return(readRDS(cache))
  env <- load_script_prefix("code/81_prior_year_pollution_cox_svi.R", "result_paths")
  dat <- as.data.table(env$analysis_dat)
  saveRDS(dat, cache, compress = FALSE)
  dat
}

spell_ids <- function(start, end) {
  stopifnot(length(start) == length(end), all(end >= start), !is.unsorted(start))
  if (!length(start)) return(integer())
  prior_end <- c(-Inf, head(cummax(as.numeric(end)), -1L))
  as.integer(cumsum(as.numeric(start) > prior_end))
}

# A listing ending on a date remains present immediately before that day's event.
# New listings affect subsequent (start, stop] risk intervals, never prior time.
concurrent_intervals <- function(start, end, index, terminal) {
  start <- as.numeric(start); end <- as.numeric(end)
  index <- as.numeric(index); terminal <- as.numeric(terminal)
  stopifnot(terminal >= index, all(end >= start))
  keep <- end > start & end > index & start < terminal
  start <- start[keep]; end <- end[keep]
  boundaries <- sort(unique(c(index, terminal, start[start > index & start < terminal],
                              end[end > index & end < terminal])))
  if (length(boundaries) < 2L) return(data.table())
  left <- head(boundaries, -1L); right <- tail(boundaries, -1L)
  count <- findInterval(left, sort(start)) - findInterval(left, sort(end))
  data.table(tstart = left - index, tstop = right - index,
             concurrent_listings = count, additional_listings = pmax(count - 1L, 0L),
             no_open_listing = as.integer(count == 0L))
}

risk_counts <- function(followup, times) {
  stopifnot(all(is.finite(followup)), all(followup >= 0))
  vapply(times, function(t) sum(followup >= t), integer(1))
}

pollutants <- c(pm25 = "pm25_prior_5ug", no2 = "no2_prior_10ppb", o3 = "o3_prior_10ppb")
organ_names <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")
baseline_adjustment <- c("age", "sex", "race", "zcta_svi_proxy", "strata(listing_center)")
eligible <- function(dat, org, term) {
  needed <- c("followup_days", "adverse_event", term, "age", "sex", "race", "zcta_svi_proxy", "listing_center")
  d <- as.data.table(dat)[WL_ORG == org & followup_days > 0]
  d <- d[complete.cases(d[, ..needed]) & is.finite(get(term))]
  droplevels(as.data.frame(d))
}
effect_row <- function(fit, term, analysis, org, pollutant, dat) {
  z <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  z <- z[z$term == term, ]
  stopifnot(nrow(z) == 1L)
  data.frame(analysis, organ = org, pollutant, n = nrow(dat),
    people = length(unique(dat$PERS_ID)), events = sum(dat$adverse_event),
    hazard_ratio = z$estimate, conf_low = z$conf.low, conf_high = z$conf.high,
    p_value = z$p.value)
}
