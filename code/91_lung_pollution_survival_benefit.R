#!/usr/bin/env Rscript
# Run from the repository root. Restricted individual data stay in ignored cache/.
source("code/r_runtime.R")
ensure_user_library()
suppressPackageStartupMessages({
  library(data.table); library(haven); library(arrow); library(dplyr)
  library(survival); library(ggplot2)
})
source("code/saf_paths.R")
source("code/daily_pollution_inputs.R")
source("code/lung_survival_benefit_helpers.R")
setDTthreads(2)

dest <- Sys.getenv("LUNG_BENEFIT_OUTPUT", "output/lung_pollution_survival_benefit")
cache <- file.path(dest, "cache")
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) message(format(Sys.time(), "%F %T"), " | ", paste0(..., collapse = ""))
bootstrap_n <- as.integer(Sys.getenv("LUNG_BENEFIT_BOOTSTRAPS", "0"))
seed <- as.integer(Sys.getenv("LUNG_BENEFIT_SEED", "20261102"))
horizon <- as.integer(Sys.getenv("LUNG_BENEFIT_HORIZON_DAYS", "365"))
stopifnot(is.finite(bootstrap_n), bootstrap_n >= 0L, is.finite(horizon), horizon > 0)
entry_start <- as.Date("2006-01-01")
las_end <- as.Date("2023-03-08")
post_end <- as.Date("2024-12-31")
prep_path <- "output/revision_2026_checklist/cache/las_era_lung/prepared.rds"
if (!file.exists(prep_path)) stop("First run code/00_prepare_lung_primary_inputs.R to prepare the established LAS-era cohort inputs")
if (identical(Sys.getenv("LUNG_BENEFIT_USE_PREPARED", "0"), "1")) {
  log_msg("Explicitly using existing prepared.rds; source data will not be refreshed")
  prepared <- readRDS(file.path(cache, "prepared.rds"))
  pre_d <- prepared$pre; post_d <- prepared$post; decision <- prepared$reference
  rm(prepared); gc(FALSE)
} else {
paths <- get_saf_paths(release = "q1_2026")
source_files <- c(prep_path, file.path(paths$pubsaf_dir, c("cand_thor.sas7bdat", "txf_lu.sas7bdat")))
stopifnot(all(file.exists(source_files)))
fwrite(data.table(file = basename(source_files), bytes = file.info(source_files)$size,
  modified = as.character(file.info(source_files)$mtime)), file.path(dest, "source_inventory.csv"))

log_msg("Reading LAS-era cohort and lung mortality/transplant sources")
prep <- readRDS(prep_path)
b <- copy(prep$cohort)
flow <- list(data.table(step = "Existing deduplicated LAS-era lung cohort", people = nrow(b)))
keep_step <- function(label, condition) {
  b <<- b[condition]
  flow[[length(flow) + 1L]] <<- data.table(step = label, people = nrow(b))
}
keep_step("Adult listings from 2006 through LAS transition", b$age >= 18 & b$index_date >= entry_start & b$index_date <= las_end)
keep_step("No other-organ candidacy in prior cohort flag", b$multi_organ_candidate == 0L)
b <- b[, .(PERS_ID, PX_ID, candidate_zip, listing_date = index_date, listing_age = age,
  sex = as.character(sex), race = as.character(race), zcta_svi_proxy,
  listing_center = as.character(listing_center), prior_waitlist_end = original_end_date)]
b[, candidate_zip := sprintf("%05d", as.integer(candidate_zip))]
cand <- as.data.table(read_sas(source_files[2], col_select = tidyselect::all_of(c(
  "PERS_ID", "PX_ID", "WL_ORG", "REC_TX_DT", "CAN_PREV_LU", "CAN_PREV_HL",
  "PERS_SSA_DEATH_DT", "PERS_OPTN_DEATH_DT", "CAN_DEATH_DT", "CAN_REM_CD", "CAN_REM_DT"))))
cand <- cand[WL_ORG == "LU"]
txf <- as.data.table(read_sas(source_files[3], col_select = tidyselect::all_of(c(
  "PERS_ID", "PX_ID", "TRR_ID", "REC_TX_DT", "TFL_PX_STAT", "TFL_PX_STAT_DT"))))
prior <- cand[, .(PX_ID, previous_lung = CAN_PREV_LU > 0 | CAN_PREV_HL > 0)]
stopifnot(!anyDuplicated(prior$PX_ID), !anyDuplicated(b$PERS_ID))
b <- merge(b, prior, by = "PX_ID", all.x = TRUE)
b <- b[!is.na(previous_lung) & previous_lung == FALSE]
flow[[length(flow) + 1L]] <- data.table(step = "No prior lung/heart-lung transplant at index registration", people = nrow(b))
# Earliest documented lung transplant, including a receiving-center registration.
# Do not require the first listing's PX_ID to equal the receiving registration.
tx_dates <- rbindlist(list(cand[!is.na(REC_TX_DT), .(PERS_ID, REC_TX_DT)],
                          txf[!is.na(REC_TX_DT), .(PERS_ID, REC_TX_DT)]))[
  , .(transplant_date = min(REC_TX_DT)), by = PERS_ID]
transplant_exits <- cand[CAN_REM_CD %in% c(4L, 14L, 15L, 18L, 19L, 21L, 22L, 23L) & !is.na(CAN_REM_DT),
  .(transplant_exit_date = min(CAN_REM_DT)), by = PERS_ID]
removal_death <- cand[CAN_REM_CD == 8L & !is.na(CAN_REM_DT), .(PERS_ID, death = CAN_REM_DT)]
deaths <- rbindlist(c(lapply(c("PERS_OPTN_DEATH_DT", "PERS_SSA_DEATH_DT", "CAN_DEATH_DT"), function(v)
  cand[!is.na(get(v)), .(PERS_ID, death = get(v))]), list(removal_death,
  txf[TFL_PX_STAT == "D" & !is.na(TFL_PX_STAT_DT), .(PERS_ID, death = TFL_PX_STAT_DT)])))
death_audit <- deaths[, .(death_date = min(death), source_dates = uniqueN(death),
  death_date_spread_days = as.numeric(max(death) - min(death))), by = PERS_ID]
alive <- txf[TFL_PX_STAT %in% c("A", "R") & !is.na(TFL_PX_STAT_DT),
  .(last_alive_date = max(TFL_PX_STAT_DT)), by = PERS_ID]
b <- Reduce(function(x, y) merge(x, y, by = "PERS_ID", all.x = TRUE), list(b, tx_dates, transplant_exits, death_audit, alive))
b <- b[(is.na(transplant_date) | transplant_date >= listing_date) &
       (is.na(death_date) | death_date >= listing_date)]
flow[[length(flow) + 1L]] <- data.table(step = "Exclude documented transplant/death before listing", people = nrow(b))
invalid_tx <- b[!is.na(transplant_date) & !is.na(death_date) & death_date < transplant_date, .N]
b[!is.na(death_date) & !is.na(transplant_date) & death_date < transplant_date, transplant_date := as.Date(NA)]
# Restrict the receiving cohort to LAS-era transplantation; later transplantation
# still censors pretransplant follow-up if it occurs before its administrative end.
b[, waiting_days := as.numeric(transplant_date - listing_date)]
pre <- copy(b)
pre[, `:=`(index_date = listing_date, age = listing_age)]
pre[, undated_transplant_censor := fifelse(is.na(transplant_date), transplant_exit_date, as.Date(NA))]
pre[, end_date := lung_date_min(death_date, transplant_date, undated_transplant_censor, rep(las_end, .N))]
# A same-day transplant/death belongs to posttransplant mortality, never both.
pre[, death_event := as.integer(!is.na(death_date) & death_date <= las_end &
  (is.na(transplant_date) | death_date < transplant_date) &
  (is.na(undated_transplant_censor) | death_date < undated_transplant_censor))]
pre[, phase := "pretransplant"]
post <- copy(b[!is.na(transplant_date) & transplant_date <= las_end])
post[, transplant_age := listing_age + waiting_days / 365.25]
post[, `:=`(index_date = transplant_date, age = transplant_age)]
# Death linkage can establish mortality after the last follow-up form. Otherwise
# censor at last known alive; never assume every recipient lived to extraction.
post[, end_date := fifelse(!is.na(death_date), pmin(death_date, post_end),
  pmin(last_alive_date, post_end))]
post[, death_event := as.integer(!is.na(death_date) & death_date <= post_end)]
no_post_followup <- post[is.na(end_date) | end_date < index_date, .N]
post <- post[!is.na(end_date) & end_date >= index_date]
post[, phase := "posttransplant"]
flow[[length(flow) + 1L]] <- data.table(step = "LAS-era recipients with usable mortality/follow-up", people = nrow(post))
fwrite(rbindlist(flow), file.path(dest, "cohort_flow.csv"))
fwrite(data.table(measure = c("People with discordant death dates", "Death before documented transplant",
  "Recipients without usable posttransplant follow-up", "Pretransplant deaths after prior waitlist exit",
  "Same-day transplant deaths", "Recipients whose last alive date follows recorded death",
  "Transplant-coded exits without linked transplant date (pretransplant censoring only)"),
  value = c(b[source_dates > 1L, .N], invalid_tx, no_post_followup,
    pre[death_event == 1L & death_date > prior_waitlist_end, .N],
    post[death_event == 1L & death_date == transplant_date, .N],
    post[!is.na(death_date) & !is.na(last_alive_date) & last_alive_date > death_date, .N],
    pre[!is.na(undated_transplant_censor), .N])),
  file.path(dest, "linkage_audit.csv"))

log_msg("Building nonoverlapping monthly follow-up and attaching official LAS")
pre_d <- lung_month_intervals(pre)
pre_d <- lung_attach_las(pre_d, prep$scores, prep$imputation)
decision <- copy(post)
decision[, start_date := transplant_date]
decision <- lung_attach_las(decision, prep$scores, prep$imputation)
post[, c("las_per10", "las_unavailable", "las_observed", "las_source_date") :=
  decision[, .(las_per10, las_unavailable, las_observed, las_source_date)]]
post_d <- lung_month_intervals(post)

log_msg("Reading monthly pollution surfaces (existing content-keyed caches where available)")
pm <- as.data.table(read_daily_pollution_aggregate("data/release/lghap_pm25_zcta_daily_parquet",
  "pm25_ug_m3", "pm25", "monthly", "output/timevarying_pollution_cox_svi/cache/pm25_daily_derived_monthly_zcta.csv.gz"))
oz <- as.data.table(read_daily_pollution_aggregate("data/release/o3_zcta_daily_parquet",
  "o3_ppb", "o3", "monthly", "output/timevarying_pollution_cox_svi/cache/o3_daily_derived_monthly_zcta.csv.gz"))
nf <- list.files("data/release/no2_zcta_monthly_parquet", "[.]parquet$", full.names = TRUE)
if (!length(nf)) stop("Monthly NO2 release is missing")
nm <- as.data.table(collect(select(open_dataset(nf), zip, year, month, no2_ppbv)))
setnames(nm, "no2_ppbv", "no2")
annual_path <- "data/release/air_pollution_zcta_parquet/air_pollution_zcta_no2_annual_2005_2025.parquet"
na <- as.data.table(read_parquet(annual_path, col_select = c("zip", "year", "no2")))
for (z in list(pm, oz, nm, na)) z[, zip := sprintf("%05d", as.integer(zip))]
monthly <- merge(pm, oz, by = c("zip", "year", "month"), all = TRUE)
monthly <- merge(monthly, nm, by = c("zip", "year", "month"), all = TRUE)
first_no2 <- min(nm$year * 12L + nm$month - 1L)
# Query grid comes from PM/O3 and monthly NO2, covering all study months.
monthly[, annual_year := year + as.integer(month == 12L) - 1L]
setnames(na, c("year", "no2"), c("annual_year", "annual_no2"))
monthly <- merge(monthly, na, by = c("zip", "annual_year"), all.x = TRUE)
monthly[, no2_annual_fallback := year * 12L + month - 1L < first_no2]
monthly[no2_annual_fallback == TRUE, no2 := annual_no2]
monthly[, exposure_source_end := as.Date(sprintf("%04d-%02d-01",
  (year * 12L + month) %/% 12L, (year * 12L + month) %% 12L + 1L)) - 1]
stopifnot(!anyDuplicated(monthly[, .(zip, year, month)]))
pre_d <- lung_attach_exposure(pre_d, monthly)
post_d <- lung_attach_exposure(post_d, monthly)
decision <- lung_attach_exposure(decision, monthly)
decision[, transplant_age := listing_age + waiting_days / 365.25]
rm(cand, txf, deaths, death_audit, alive, monthly, pm, oz, nm, na, prep); gc(FALSE)
saveRDS(list(pre = pre_d, post = post_d, reference = decision), file.path(cache, "prepared.rds"), compress = FALSE)
}

pollutants <- c("pm25", "no2", "o3")
models <- list(); result <- list(); model_log <- list(); curves <- list(); refs <- list(); grids <- list()
for (p in pollutants) {
  models[[p]] <- list()
  for (phase_name in c("pretransplant", "posttransplant")) {
    d <- if (phase_name == "pretransplant") pre_d else post_d
    log_msg("Fitting ", phase_name, " / ", p)
    m <- lung_fit_cox(d, p)
    models[[p]][[phase_name]] <- m
    s <- summary(m$fit)$coefficients["exposure_scaled", ]
    ci <- confint(m$fit)["exposure_scaled", ]
    result[[length(result) + 1L]] <- data.table(phase = phase_name, pollutant = p,
      people = uniqueN(m$data$PERS_ID), intervals = nrow(m$data), deaths = sum(m$data$event),
      person_years = sum(m$data$tstop - m$data$tstart) / 365.25,
      HR = exp(s["coef"]), HR_low = exp(ci[1]), HR_high = exp(ci[2]), p_value = s["Pr(>|z|)"],
      no2_fallback_person_years = if (p == "no2") sum((m$data$tstop - m$data$tstart)[m$data$no2_annual_fallback]) / 365.25 else 0)
    model_log[[length(model_log) + 1L]] <- data.table(phase = phase_name, pollutant = p,
      formula = m$formula, status = if (length(m$warnings)) paste(m$warnings, collapse = "; ") else "OK")
    saveRDS(m$fit, file.path(cache, paste0(phase_name, "_", p, "_cox.rds")))
  }
  a <- models[[p]]$pretransplant; z <- models[[p]]$posttransplant
  r <- decision[PERS_ID %in% a$data$PERS_ID & PERS_ID %in% z$data$PERS_ID &
    complete.cases(decision[, .(listing_age, transplant_age, sex, race, zcta_svi_proxy, listing_center, las_per10, las_unavailable)]) &
    is.finite(get(p))]
  r <- lung_supported_reference(r, a, z, horizon)
  if (!nrow(r)) stop("No supported paired recipients for ", p)
  refs[[p]] <- r
  lo <- max(quantile(a$data[[p]], .1), quantile(z$data[[p]], .1))
  hi <- min(quantile(a$data[[p]], .9), quantile(z$data[[p]], .9))
  if (lo >= hi) stop("No common central exposure range for ", p)
  grids[[p]] <- seq(lo, hi, length.out = 21L)
  curves[[p]] <- lung_predict_benefit(a, z, r, p, grids[[p]], horizon)
}
fwrite(rbindlist(result), file.path(dest, "pollutant_results.csv"))
fwrite(rbindlist(model_log), file.path(dest, "model_log.csv"))
point <- rbindlist(curves)
fwrite(point, file.path(dest, "survival_benefit_curves.csv"))
saveRDS(list(models = models, reference = refs, grids = grids, point = point),
  file.path(cache, "fitted.rds"), compress = FALSE)

# Everything below can also be rerun against the local fitted cache with
# code/92_lung_pollution_benefit_bootstrap.R, without rereading restricted SAFs.
# Resample PEOPLE jointly across both phases. Multiplicity weights avoid physically
# expanding each resample. Robust SEs are unnecessary inside the bootstrap.
set.seed(seed)
ids <- sort(unique(c(pre_d$PERS_ID, post_d$PERS_ID)))
boot <- list(); boot_log <- list()
# Replace old uncertainty outputs even on a zero-bootstrap/debug rerun.
fwrite(data.table(iteration = integer(), pollutant = character(), status = character()), file.path(dest, "bootstrap_log.csv"))
empty_boot <- copy(point[0]); empty_boot[, iteration := integer()]
fwrite(empty_boot, file.path(dest, "bootstrap_curves.csv"))
if (bootstrap_n > 0L) for (iteration in seq_len(bootstrap_n)) {
  counts <- tabulate(sample.int(length(ids), length(ids), replace = TRUE), nbins = length(ids))
  weights <- setNames(counts, as.character(ids))
  for (p in pollutants) {
    log_msg("Bootstrap ", iteration, "/", bootstrap_n, " / ", p)
    ans <- tryCatch({
      a <- lung_fit_cox(models[[p]]$pretransplant$data, p, weights)
      z <- lung_fit_cox(models[[p]]$posttransplant$data, p, weights)
      # Fixed empirical recipient population: uncertainty in fitted survival models.
      # Require every original prediction stratum to remain represented.
      needed <- unique(as.character(refs[[p]]$listing_center))
      if (!all(needed %in% a$data$listing_center) || !all(needed %in% z$data$listing_center)) stop("Bootstrap omitted reference center")
      q <- lung_predict_benefit(a, z, refs[[p]], p, grids[[p]], horizon)
      q[, iteration := iteration]
      q
    }, error = function(e) e)
    status <- if (inherits(ans, "error")) conditionMessage(ans) else "OK"
    boot_log[[length(boot_log) + 1L]] <- data.table(iteration, pollutant = p, status)
    if (!inherits(ans, "error")) boot[[length(boot) + 1L]] <- ans
  }
  fwrite(rbindlist(boot_log), file.path(dest, "bootstrap_log.csv"))
  if (length(boot)) fwrite(rbindlist(boot), file.path(dest, "bootstrap_curves.csv"))
}
metrics <- c("survival_without_transplant", "survival_after_transplant", "benefit_probability")
if (length(boot)) {
  bd <- rbindlist(boot)
  bands <- bd[, c(list(successful_bootstraps = uniqueN(iteration)),
    unlist(lapply(.SD, function(x) list(low = quantile(x, .025), high = quantile(x, .975))), recursive = FALSE)),
    by = .(pollutant, concentration), .SDcols = metrics]
  setnames(bands, old = setdiff(names(bands), c("pollutant", "concentration", "successful_bootstraps")),
    new = as.vector(t(outer(metrics, c("low", "high"), paste, sep = "_"))))
  point <- merge(point, bands, by = c("pollutant", "concentration"), all.x = TRUE)
  point[is.na(successful_bootstraps), successful_bootstraps := 0L]
} else point[, successful_bootstraps := 0L]
point[, confidence_status := fifelse(successful_bootstraps == 0L, "not_run_or_all_failed",
  fifelse(successful_bootstraps >= max(100L, ceiling(.8 * bootstrap_n)), "estimated", "insufficient_resamples"))]
ci_cols <- grep("_(low|high)$", names(point), value = TRUE)
for (nm in ci_cols) point[confidence_status != "estimated", (nm) := NA_real_]
fwrite(point, file.path(dest, "survival_benefit_curves.csv"))
contrast <- rbindlist(lapply(pollutants, function(p) {
  q <- point[pollutant == p][order(concentration)]
  data.table(pollutant = p, low_exposure = q$concentration[1], high_exposure = tail(q$concentration, 1),
    benefit_difference_pp = 100 * (tail(q$benefit_probability, 1) - q$benefit_probability[1]),
    recipients = q$recipients[1])
}))
if (length(boot)) {
  bc <- rbindlist(boot)[order(concentration), .(delta = 100 * (tail(benefit_probability, 1) - benefit_probability[1])), by = .(pollutant, iteration)]
  contrast <- merge(contrast, bc[, .(low_pp = quantile(delta, .025), high_pp = quantile(delta, .975),
    successful_bootstraps = .N), by = pollutant], by = "pollutant", all.x = TRUE)
}
contrast <- merge(contrast, unique(point[, .(pollutant, confidence_status)]), by = "pollutant")
if ("low_pp" %in% names(contrast)) contrast[confidence_status != "estimated", c("low_pp", "high_pp") := .(NA_real_, NA_real_)]
fwrite(contrast, file.path(dest, "benefit_exposure_contrasts.csv"))
fig <- ggplot(point, aes(concentration, 100 * benefit_probability)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_line(linewidth = .9, color = "#226078") +
  facet_wrap(~pollutant, scales = "free_x", labeller = as_labeller(c(pm25 = "PM2.5 (ug/m3)", no2 = "NO2 (ppb)", o3 = "O3 (ppb)"))) +
  labs(x = "Previous complete month's concentration", y = "Estimated transplant survival benefit (percentage points)",
    title = paste0("Model-based lung transplant benefit at ", horizon, " days"),
    subtitle = "Posttransplant survival minus conditional survival without transplant",
    caption = paste("Standardized to supported LAS-era recipients; covariates frozen at transplantation. Observational model estimates.",
      if (all(point$confidence_status == "estimated")) "Shading: pointwise 95% bootstrap intervals." else "Point estimates; bootstrap uncertainty incomplete.")) +
  theme_bw(base_size = 11)
if ("benefit_probability_low" %in% names(point)) fig <- fig +
  geom_ribbon(aes(ymin = 100 * benefit_probability_low, ymax = 100 * benefit_probability_high),
    fill = "#226078", alpha = .15)
ggsave(file.path(dest, "survival_benefit.png"), fig, width = 12, height = 5, dpi = 170)
ggsave(file.path(dest, "survival_benefit.pdf"), fig, width = 12, height = 5)
writeLines(c(
  "Two separate time-varying cause-specific Cox mortality models; no target-trial emulation or treatment/censoring weights.",
  "Adults, first documented lung transplant candidates, no other-organ candidacy per existing cohort flag; listing entry 2006-01-01 through 2023-03-08.",
  "Pretransplant follow-up ends at transplant, death or 2023-03-08, including deaths after delisting. Negative death ascertainment assumes registry linkage through this date.",
  "Posttransplant follow-up ends at documented death, last known alive, or 2024-12-31; all-cause survival continues through retransplantation when follow-up is available.",
  "Earliest available death date across OPTN, SSA, candidate death/removal and deceased follow-up; discordant dates audited. Same-day transplant death contributes 0.5 postoperative day.",
  "Transplant-coded removal without any linked transplant date censors pretransplant mortality follow-up at removal; it cannot initiate posttransplant follow-up.",
  "Prior complete calendar-month exposure updates monthly. Before monthly NO2 starts, use previous DECISION-YEAR annual NO2; never same-year annual fallback. No filling gaps in the monthly era.",
  "Permanent listing ZIP and listing-year ACS ZCTA vulnerability are retained in both phases; residential relocation is not captured.",
  "Same prior lung adjustment: age, sex, race, zcta_svi_proxy, official LAS per 10 with missing indicator, strata(listing_center). Constant covariates omitted. No new calendar/donor adjustment.",
  "Pretransplant age is age at listing, LAS updates at monthly starts. Posttransplant age is age at transplant and LAS is frozen at transplant; missing LAS uses prior study's reference median.",
  "Efron ties, person-clustered robust coefficient intervals. PM2.5 HR per 5 ug/m3; NO2/O3 per 10 ppb. Deterioration alone is not mortality.",
  "Pretransplant survival prediction: exp(-(H0(waiting_days+horizon)-H0(waiting_days))*exp(lp)); posttransplant: exp(-H0(horizon)*exp(lp)).",
  "Predictions hold exposure and LAS fixed, preserve each phase's age definition and listing-center stratum, and average both predictions over identical recipients within each pollutant.",
  "Prediction exposure grid spans the overlap of the two models' central 10th-90th percentile interval distributions; it is not a regulatory threshold.",
  "Require at least 5 observed people per center still under follow-up at each prediction horizon in each phase. This is a support check, not causal identification.",
  paste0("Horizon days: ", horizon, "; requested joint person bootstrap resamples: ", bootstrap_n, "; seed: ", seed, ". Percentile intervals are pointwise and conditional on the fixed recipient reference population."),
  "Check bootstrap_log.csv for failures; zero or few successful replicates do not provide reliable uncertainty. Script 91 defaults to point estimates (0 resamples); script 92 defaults to 200 for uncertainty estimation.",
  "Benefit is an adjusted model-based survival difference, not an identified causal treatment or pollution effect. No allocation policy changes are evaluated."
), file.path(dest, "methods.txt"))
capture.output(sessionInfo(), file = file.path(dest, "session_info.txt"))
log_msg("Completed: ", dest)
print(rbindlist(result)); print(contrast)
