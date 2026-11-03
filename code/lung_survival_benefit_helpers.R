# Helpers for descriptive, model-based lung transplant survival benefit.
# Source after loading data.table and survival. No source-data reads or writes here.

lung_date_min <- function(...) {
  z <- do.call(pmin, c(lapply(list(...), as.numeric), na.rm = TRUE))
  z[!is.finite(z)] <- NA_real_
  as.Date(z, origin = "1970-01-01")
}

lung_month_intervals <- function(cohort) {
  d <- copy(as.data.table(cohort))
  d <- d[!is.na(index_date) & !is.na(end_date) & end_date >= index_date]
  if (!nrow(d)) return(d)
  # Day-zero deaths contribute half a day; ordinary follow-up is (start, stop].
  # Calendar-month boundaries are shared, never duplicated or overlapped.
  first <- as.integer(format(d$index_date, "%Y")) * 12L + as.integer(format(d$index_date, "%m")) - 1L
  last <- as.integer(format(d$end_date, "%Y")) * 12L + as.integer(format(d$end_date, "%m")) - 1L
  lengths <- last - first + 1L
  ids <- rep(first, lengths) + sequence(lengths) - 1L
  out <- d[rep(seq_len(nrow(d)), lengths)]
  month_date <- function(x) as.Date(sprintf("%04d-%02d-01", x %/% 12L, x %% 12L + 1L))
  out[, start_date := pmax(index_date, month_date(ids))]
  out[, stop_date := pmin(end_date, month_date(ids + 1L))]
  out[, `:=`(tstart = as.numeric(start_date - index_date),
             tstop = as.numeric(stop_date - index_date),
             event = as.integer(death_event == 1L & stop_date == end_date))]
  # A death on a new month's first day belongs to the preceding interval.
  out <- out[tstop > tstart | (tstart == 0 & tstop == 0 & event == 1L)]
  out[tstop == 0 & event == 1L, tstop := 0.5]
  setorder(out, PERS_ID, tstart)
  out[]
}

lung_attach_exposure <- function(d, monthly) {
  out <- copy(d)
  # Previous complete month. Before monthly NO2 begins, the reader supplies
  # the PREVIOUS calendar year's annual value (never a future-containing year).
  id <- as.integer(format(out$start_date, "%Y")) * 12L + as.integer(format(out$start_date, "%m")) - 2L
  out[, `:=`(exposure_year = id %/% 12L, exposure_month = id %% 12L + 1L)]
  m <- copy(monthly)
  setkey(m, zip, year, month)
  cols <- c("pm25", "no2", "o3", "no2_annual_fallback", "exposure_source_end")
  values <- m[out, on = .(zip = candidate_zip, year = exposure_year, month = exposure_month), ..cols]
  for (nm in cols) set(out, j = nm, value = values[[nm]])
  stopifnot(all(out$exposure_source_end < out$start_date, na.rm = TRUE))
  out[]
}

lung_attach_las <- function(d, scores, imputation) {
  out <- copy(d)
  s <- copy(scores)
  setkey(s, PERS_ID, score_date)
  out[, las_observed := s[.SD, on = .(PERS_ID, score_date = start_date), roll = Inf, x.las]]
  out[, las_source_date := s[.SD, on = .(PERS_ID, score_date = start_date), roll = Inf, x.score_date]]
  out[, las_unavailable := as.integer(!is.finite(las_observed))]
  out[, las_per10 := fifelse(las_unavailable == 1L, imputation, las_observed) / 10]
  stopifnot(all(out$las_source_date <= out$start_date, na.rm = TRUE))
  out[]
}

lung_fit_cox <- function(d, pollutant, weights = NULL, initial = NULL) {
  cols <- c("PERS_ID", "tstart", "tstop", "event", pollutant, "age", "sex", "race",
            "zcta_svi_proxy", "listing_center", "las_per10", "las_unavailable")
  keep <- unique(c(cols, "no2_annual_fallback"))
  keep <- intersect(keep, names(d))
  x <- copy(d[complete.cases(d[, ..cols]) & is.finite(get(pollutant)) & tstop > tstart, ..keep])
  x[, exposure_scaled := get(pollutant) / if (pollutant == "pm25") 5 else 10]
  if (is.null(weights)) x[, fit_weight := 1] else {
    x[, fit_weight := unname(weights[as.character(PERS_ID)])]
    x <- x[is.finite(fit_weight) & fit_weight > 0]
  }
  if (sum(x$event) < 20L) stop("Fewer than 20 deaths for ", pollutant)
  # Match the prior model's adjustment; omit only genuinely constant columns.
  adjust <- c("age", "sex", "race", "zcta_svi_proxy", "las_per10", "las_unavailable")
  adjust <- adjust[vapply(adjust, function(v) uniqueN(x[[v]]) > 1L, logical(1))]
  # Do not serialize this function's frame (including full input datasets) in
  # every fitted formula. All variables are supplied through model data.
  form <- as.formula(paste("Surv(tstart, tstop, event) ~",
    paste(c("exposure_scaled", adjust, "strata(listing_center)"), collapse = " + ")), env = globalenv())
  columns <- colnames(model.matrix(reformulate(c("exposure_scaled", adjust)), data = as.data.frame(x)))
  columns <- setdiff(columns, "(Intercept)")
  if (is.null(initial)) initial <- rep(0, length(columns)) else {
    initial <- unname(initial[columns])
    initial[!is.finite(initial)] <- 0
  }
  warnings <- character()
  fit <- withCallingHandlers(if (is.null(weights)) {
    coxph(form, data = as.data.frame(x), ties = "efron", weights = fit_weight, init = initial,
      cluster = PERS_ID, robust = TRUE, model = TRUE, x = FALSE, y = TRUE)
  } else {
    coxph(form, data = as.data.frame(x), ties = "efron", weights = fit_weight, init = initial,
      robust = FALSE, model = TRUE, x = FALSE, y = TRUE)
  }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  if (any(!is.finite(coef(fit)))) stop("Nonfinite coefficient for ", pollutant)
  list(fit = fit, data = x, warnings = unique(warnings), formula = paste(deparse(form), collapse = " "))
}

lung_baseline <- function(model) {
  h <- as.data.table(basehaz(model$fit, centered = FALSE))
  if (!"strata" %in% names(h)) h[, strata := unique(as.character(model$data$listing_center))]
  h[, strata := as.character(strata)]
  h
}

lung_hazard_at <- function(h, centers, times) {
  ans <- rep(NA_real_, length(times))
  for (ctr in unique(centers)) {
    ids <- which(centers == ctr)
    z <- h[strata == ctr]
    if (!nrow(z)) next
    ans[ids] <- c(0, z$hazard)[findInterval(times[ids], z$time) + 1L]
  }
  ans
}

lung_supported_reference <- function(reference, pre, post, horizon = 365, min_at_risk = 5L) {
  r <- copy(reference)
  r[, supported := FALSE]
  for (ctr in unique(as.character(r$listing_center))) {
    ids <- which(as.character(r$listing_center) == ctr)
    a <- pre$data[as.character(listing_center) == ctr, .(last = max(tstop)), by = PERS_ID]$last
    b <- post$data[as.character(listing_center) == ctr, .(last = max(tstop)), by = PERS_ID]$last
    if (length(a) < min_at_risk || length(b) < min_at_risk) next
    max_pre <- sort(a, decreasing = TRUE)[min_at_risk]
    max_post <- sort(b, decreasing = TRUE)[min_at_risk]
    r[ids, supported := waiting_days + horizon <= max_pre & horizon <= max_post]
  }
  r[supported == TRUE]
}

lung_predict_benefit <- function(pre, post, reference, pollutant, grid, horizon = 365) {
  if (!nrow(reference)) stop("No recipients have common center/horizon support")
  pre_h <- lung_baseline(pre); post_h <- lung_baseline(post)
  ctr <- as.character(reference$listing_center)
  pre_delta <- lung_hazard_at(pre_h, ctr, reference$waiting_days + horizon) -
    lung_hazard_at(pre_h, ctr, reference$waiting_days)
  post_delta <- lung_hazard_at(post_h, ctr, rep(horizon, nrow(reference)))
  rbindlist(lapply(grid, function(p) {
    a <- as.data.frame(reference); b <- a
    a$age <- a$listing_age; b$age <- b$transplant_age
    a$exposure_scaled <- b$exposure_scaled <- p / if (pollutant == "pm25") 5 else 10
    # Explicit zero reference matches basehaz(centered=FALSE).
    lp_a <- predict(pre$fit, newdata = a, type = "lp", reference = "zero")
    lp_b <- predict(post$fit, newdata = b, type = "lp", reference = "zero")
    s_a <- exp(-pre_delta * exp(lp_a)); s_b <- exp(-post_delta * exp(lp_b))
    if (any(!is.finite(s_a)) || any(!is.finite(s_b))) stop("Invalid survival predictions")
    data.table(pollutant, concentration = p, horizon_days = horizon, recipients = nrow(reference),
      survival_without_transplant = mean(s_a), survival_after_transplant = mean(s_b),
      benefit_probability = mean(s_b - s_a))
  }))
}

# Exact area under the fitted Cox step-survival curve. Each candidate's waiting
# time shifts the pretransplant event knots; no trapezoidal interpolation or
# multiplication of endpoint survival by 365 is used. Return population means.
lung_rmst_integrate <- function(h, centers, starts, relative_risk, horizons) {
  relative_risk <- as.matrix(relative_risk)
  stopifnot(nrow(relative_risk) == length(starts), length(centers) == length(starts),
    length(starts) > 0L, all(is.finite(relative_risk)), all(relative_risk >= 0),
    all(is.finite(starts)), all(starts >= 0), all(is.finite(horizons)), all(horizons >= 0))
  total <- matrix(0, nrow = length(horizons), ncol = ncol(relative_risk))
  for (ctr in unique(centers)) {
    z <- h[strata == ctr]
    if (!nrow(z)) stop("Missing baseline hazard for center ", ctr)
    at <- function(t) c(0, z$hazard)[findInterval(t, z$time) + 1L]
    ids <- which(centers == ctr)
    # Date-derived Cox times lie on an integer/half-day lattice. Integrate on
    # its shared knots in matrix batches instead of rebuilding knots for every
    # waiting time. Extra knots do not approximate or alter the step-function
    # area. Arbitrary continuous-time inputs use the general exact path below.
    fractions <- unique(as.vector(outer(unique(z$time %% 1), unique(starts[ids] %% 1), "-")) %% 1)
    if (length(fractions) <= 4L) {
      knots <- sort(unique(c(0, max(horizons), as.vector(outer(seq(0, floor(max(horizons))), fractions, "+")))))
      knots <- knots[knots >= 0 & knots <= max(horizons)]
      if (length(knots) < 2L) next
      left <- head(knots, -1L); right <- tail(knots, -1L)
      widths <- outer(right, horizons, pmin) - outer(left, horizons, pmin)
      for (chunk in split(ids, ceiling(seq_along(ids) / 256L))) {
        hz <- matrix(at(as.vector(outer(left, starts[chunk], "+"))), nrow = length(left))
        delta <- pmax(sweep(hz, 2L, at(starts[chunk]), "-"), 0)
        for (j in seq_len(ncol(relative_risk))) {
          area <- crossprod(widths, exp(-sweep(delta, 2L, relative_risk[chunk, j], "*")))
          total[, j] <- total[, j] + rowSums(area)
        }
      }
      next
    }
    for (group in split(ids, starts[ids])) {
      start <- starts[group[1L]]
      knots <- sort(unique(c(0, z$time[z$time > start & z$time < start + max(horizons)] - start,
        max(horizons))))
      if (length(knots) < 2L) next
      left <- head(knots, -1L); right <- tail(knots, -1L)
      delta <- pmax(at(start + left) - at(start), 0)
      widths <- outer(right, horizons, pmin) - outer(left, horizons, pmin)
      for (chunk in split(group, ceiling(seq_along(group) / 128L))) {
        rr <- as.vector(relative_risk[chunk, , drop = FALSE])
        area <- crossprod(widths, exp(-outer(delta, rr)))
        for (j in seq_len(ncol(relative_risk))) {
          cols <- (j - 1L) * length(chunk) + seq_along(chunk)
          total[, j] <- total[, j] + rowSums(area[, cols, drop = FALSE])
        }
      }
    }
  }
  t(total / length(starts))
}

lung_predict_rmst <- function(pre, post, reference, pollutant, grid, horizons = 365,
                              baselines = NULL) {
  if (!nrow(reference)) stop("No reference recipients for RMST")
  if (is.null(baselines)) baselines <- list(pre = lung_baseline(pre), post = lung_baseline(post))
  a <- b <- as.data.frame(reference)
  a$age <- a$listing_age; b$age <- b$transplant_age
  # Exposure is additive in the existing models; construct all scenarios from
  # the zero-exposure linear predictor while preserving phase-specific age.
  a$exposure_scaled <- b$exposure_scaled <- 0
  scale <- if (pollutant == "pm25") 5 else 10
  lp_a <- predict(pre$fit, newdata = a, type = "lp", reference = "zero")
  lp_b <- predict(post$fit, newdata = b, type = "lp", reference = "zero")
  rr_a <- exp(outer(lp_a, grid / scale * coef(pre$fit)["exposure_scaled"], "+"))
  rr_b <- exp(outer(lp_b, grid / scale * coef(post$fit)["exposure_scaled"], "+"))
  centers <- as.character(reference$listing_center)
  without <- lung_rmst_integrate(baselines$pre, centers, reference$waiting_days, rr_a, horizons)
  after <- lung_rmst_integrate(baselines$post, centers, rep(0, nrow(reference)), rr_b, horizons)
  out <- data.table(pollutant, concentration = rep(grid, times = length(horizons)),
    horizon_days = rep(horizons, each = length(grid)), recipients = nrow(reference),
    rmst_without_transplant_days = as.vector(without), rmst_after_transplant_days = as.vector(after),
    rmst_benefit_days = as.vector(after - without))
  stopifnot(all(out$rmst_without_transplant_days >= -1e-8),
    all(out$rmst_after_transplant_days >= -1e-8),
    all(out$rmst_without_transplant_days <= out$horizon_days + 1e-8),
    all(out$rmst_after_transplant_days <= out$horizon_days + 1e-8))
  out[]
}

lung_validate_rmst <- function(curves, times, endpoint) {
  x <- merge(curves, as.data.table(endpoint), by = c("pollutant", "concentration"))
  end <- merge(times[horizon_days == 365], curves,
    by = c("pollutant", "concentration"), suffixes = c("_time", "_exposure"))
  checks <- c(
    complete_exposure_grid = nrow(x) == nrow(curves),
    rmst_exceeds_endpoint_survival_area = all(x$rmst_without_transplant_days >= 365 * x$survival_without_transplant - 1e-8) &&
      all(x$rmst_after_transplant_days >= 365 * x$survival_after_transplant - 1e-8),
    bounded_by_365_days = all(curves$rmst_without_transplant_days <= 365 + 1e-8) && all(curves$rmst_after_transplant_days <= 365 + 1e-8),
    time_and_exposure_endpoints_agree = nrow(end) == 6L &&
      all(abs(end$rmst_benefit_days_time - end$rmst_benefit_days_exposure) < 1e-8),
    zero_area_at_time_zero = all(abs(times[horizon_days == 0, rmst_benefit_days]) < 1e-8),
    increasing_rmst_over_time = all(times[order(horizon_days),
      .(ok = all(diff(rmst_without_transplant_days) >= -1e-8) && all(diff(rmst_after_transplant_days) >= -1e-8)),
      by = .(pollutant, concentration)]$ok))
  if (!all(checks)) stop("RMST output validation failed: ", paste(names(checks)[!checks], collapse = ", "))
  data.table(check = names(checks), passed = unname(checks))
}
