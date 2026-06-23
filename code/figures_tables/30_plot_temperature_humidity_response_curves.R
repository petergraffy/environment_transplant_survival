#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(haven)
  library(readr)
  library(splines)
  library(stringr)
  library(survival)
  library(tibble)
  library(tidyr)
})

saf_paths_source <- file.path("code", "saf_paths.R")
if (!file.exists(saf_paths_source)) saf_paths_source <- "saf_paths.R"
source(saf_paths_source)
saf_paths <- get_saf_paths()
assert_saf_files(saf_paths, include_stathist = TRUE)

pubsaf_dir <- saf_paths$pubsaf_dir
supp_dir <- saf_paths$supp_dir
exposure_dir <- file.path("output", "waitlist_environment_exposure")
out_dir <- file.path("output", "figures", "waitlist_environment_exposure")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analysis_start_year <- as.integer(Sys.getenv("ENV_ANALYSIS_START_YEAR", "2006"))
analysis_end_year <- as.integer(Sys.getenv("ENV_ANALYSIS_END_YEAR", "2023"))
analysis_end_date <- as.Date(sprintf("%04d-12-31", analysis_end_year))
target_organs <- c("HR", "KI", "LI", "LU")
organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

clean_zip <- function(x) {
  z <- str_extract(as.character(x), "[0-9]+")
  z <- if_else(nchar(z) > 5L, substr(z, 1, 5), z)
  if_else(is.na(z), NA_character_, str_pad(z, 5, side = "left", pad = "0"))
}

as_factor_missing <- function(x) {
  factor(if_else(is.na(x) | as.character(x) == "", "Missing", as.character(x)))
}

flag_yes <- function(x) {
  y <- str_to_upper(str_trim(as.character(x)))
  as.integer(y %in% c("1", "Y", "YES", "TRUE", "T"))
}

safe_log1p <- function(x) {
  log(pmax(as.numeric(x), 0) + 1)
}

median_impute <- function(x) {
  x <- as.numeric(x)
  med <- median(x, na.rm = TRUE)
  if (!is.finite(med)) med <- NA_real_
  if_else(is.na(x), med, x)
}

collapse_rare <- function(x, min_n = 100L) {
  y <- as.character(x)
  y[is.na(y) | y == ""] <- "Missing"
  tab <- table(y)
  y[!y %in% names(tab)[tab >= min_n]] <- "Other"
  factor(y)
}

calc_egfr_2021 <- function(creatinine, age, sex) {
  scr <- as.numeric(creatinine)
  age <- as.numeric(age)
  sex_chr <- str_to_upper(as.character(sex))
  female <- sex_chr %in% c("F", "FEMALE")
  kappa <- if_else(female, 0.7, 0.9)
  alpha <- if_else(female, -0.241, -0.302)
  sex_mult <- if_else(female, 1.012, 1.0)
  142 * pmin(scr / kappa, 1)^alpha * pmax(scr / kappa, 1)^(-1.200) * 0.9938^age * sex_mult
}

candidate_files <- tribble(
  ~candidate_group, ~path,
  "kidney_pancreas", file.path(pubsaf_dir, "cand_kipa.sas7bdat"),
  "liver_intestine", file.path(pubsaf_dir, "cand_liin.sas7bdat"),
  "thoracic", file.path(pubsaf_dir, "cand_thor.sas7bdat")
)

candidate_cols <- c(
  "PERS_ID", "PX_ID", "WL_ORG", "CAN_LISTING_DT", "CAN_ACTIVATE_DT", "CAN_SOURCE",
  "CAN_REM_DT", "CAN_REM_CD", "CAN_DEATH_DT", "CAN_ENDWLFU", "CAN_AGE_AT_LISTING",
  "CAN_GENDER", "CAN_RACE_SRTR", "CAN_LISTING_CTR_CD", "CAN_INIT_STAT",
  "CAN_ON_DIAL", "CAN_DIAL", "CAN_DIAL_DT", "CAN_MOST_RECENT_CREAT",
  "CAN_INIT_SRTR_LAB_MELD", "CAN_LAST_SRTR_LAB_MELD", "CAN_TOT_BILI",
  "CAN_TOT_ALBUMIN", "CAN_LAST_SERUM_SODIUM", "CAN_FVC", "CAN_FEV1", "CAN_PCO2",
  "CAN_AT_REST_O2", "CAN_SIX_MIN_WALK", "CAN_PULM_ART_MEAN", "CAN_CARDIAC_OUTPUT",
  "CAN_VENTILATOR", "CAN_ON_VENTILATOR", "CAN_ECMO", "CAN_VAD_TAH", "CAN_VAD1",
  "CAN_VAD2", "CAN_CORTICOST_DEPND"
)

log_msg("Reading saved waitlist-period environment exposures")
exposures <- read_csv(
  file.path(exposure_dir, "waitlist_environment_exposures.csv.gz"),
  col_types = cols(
    waitlist_row_id = col_double(),
    PERS_ID = col_double(),
    PX_ID = col_double(),
    WL_ORG = col_character(),
    candidate_zip = col_character(),
    .default = col_guess()
  ),
  show_col_types = FALSE
) %>%
  select(
    waitlist_row_id, PERS_ID, PX_ID, WL_ORG, candidate_zip,
    tmean_c_waitlist_mean, temp_range_c_waitlist_mean, rhmean_pct_waitlist_mean,
    pm25_waitlist_ug_m3, o3_waitlist_ppb, no2_waitlist
  ) %>%
  mutate(
    waitlist_row_id = as.integer(waitlist_row_id),
    candidate_zip = clean_zip(candidate_zip),
    tmean_waitlist_5c = tmean_c_waitlist_mean / 5,
    temp_range_waitlist_5c = temp_range_c_waitlist_mean / 5,
    rhmean_waitlist_10pct = rhmean_pct_waitlist_mean / 10,
    pm25_waitlist_5ug = pm25_waitlist_ug_m3 / 5,
    o3_waitlist_10ppb = o3_waitlist_ppb / 10,
    no2_waitlist_10unit = no2_waitlist / 10
  )

log_msg("Rebuilding model covariates from SAF")
candidate_data <- candidate_files %>%
  rowwise() %>%
  mutate(data = list(read_sas(path, col_select = any_of(candidate_cols)) %>% mutate(candidate_group = candidate_group))) %>%
  ungroup() %>%
  select(data) %>%
  unnest(data)

candidate_zip <- read_sas(file.path(supp_dir, "canzip2506.sas7bdat")) %>%
  transmute(PERS_ID, PX_ID, candidate_zip = clean_zip(CAN_PERM_ZIP)) %>%
  distinct(PERS_ID, PX_ID, .keep_all = TRUE)

heart_just <- read_sas(
  file.path(pubsaf_dir, "statjust_hr1a.sas7bdat"),
  col_select = any_of(c(
    "PX_ID", "WL_ORG", "CANHX_CHG_DT", "CANHX_LVAD_TYPE", "CANHX_RVAD_TYPE",
    "CANHX_ECMO", "CANHX_LAB_SERUM_CREAT", "CANHX_LAB_BILI",
    "CANHX_LAB_ALBUMIN", "CANHX_LAB_SODIUM", "CANHX_LAB_BNP"
  ))
) %>%
  filter(WL_ORG == "HR") %>%
  arrange(PX_ID, CANHX_CHG_DT) %>%
  group_by(PX_ID) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  transmute(
    PX_ID,
    hr_short_mcs = as.integer(flag_yes(CANHX_ECMO) == 1L | flag_yes(CANHX_RVAD_TYPE) == 1L),
    hr_durable_lvad = flag_yes(CANHX_LVAD_TYPE),
    hr_creatinine = CANHX_LAB_SERUM_CREAT,
    hr_bilirubin = CANHX_LAB_BILI,
    hr_albumin = CANHX_LAB_ALBUMIN,
    hr_sodium = CANHX_LAB_SODIUM,
    hr_bnp = CANHX_LAB_BNP
  )

cohort <- candidate_data %>%
  filter(WL_ORG %in% target_organs) %>%
  left_join(candidate_zip, by = c("PERS_ID", "PX_ID")) %>%
  left_join(heart_just, by = "PX_ID") %>%
  mutate(
    waitlist_row_id = row_number(),
    index_date = coalesce(CAN_ACTIVATE_DT, CAN_LISTING_DT),
    index_year = as.integer(format(index_date, "%Y")),
    event_date = coalesce(CAN_REM_DT, CAN_DEATH_DT, CAN_ENDWLFU),
    observed_end_date = pmin(event_date, analysis_end_date),
    event_observed = !is.na(event_date) & event_date <= analysis_end_date,
    adverse_event = as.integer(event_observed & CAN_REM_CD %in% c(8, 13)),
    followup_days = as.numeric(observed_end_date - index_date),
    age = CAN_AGE_AT_LISTING,
    sex = as_factor_missing(if_else(CAN_GENDER %in% c("M", "F"), CAN_GENDER, NA_character_)),
    race = collapse_rare(CAN_RACE_SRTR, min_n = 100L),
    listing_center = collapse_rare(CAN_LISTING_CTR_CD, min_n = 100L),
    index_year_centered = index_year - 2015,
    kidney_on_dialysis = flag_yes(coalesce(as.character(CAN_ON_DIAL), as.character(CAN_DIAL))),
    dialysis_years = as.numeric(index_date - CAN_DIAL_DT) / 365.25,
    dialysis_years = if_else(is.na(dialysis_years) | dialysis_years < 0, 0, dialysis_years),
    kidney_dialysis_age_score = age / 10 + 0.8 * kidney_on_dialysis + 0.3 * log1p(dialysis_years),
    liver_meld = coalesce(CAN_INIT_SRTR_LAB_MELD, CAN_LAST_SRTR_LAB_MELD),
    hr_creatinine_for_score = coalesce(hr_creatinine, CAN_MOST_RECENT_CREAT),
    hr_albumin = coalesce(hr_albumin, CAN_TOT_ALBUMIN),
    hr_bilirubin = coalesce(hr_bilirubin, CAN_TOT_BILI),
    hr_sodium = coalesce(hr_sodium, CAN_LAST_SERUM_SODIUM),
    hr_short_mcs = coalesce(hr_short_mcs, flag_yes(CAN_ECMO), 0L),
    hr_durable_lvad = coalesce(hr_durable_lvad, flag_yes(CAN_VAD_TAH), 0L)
  ) %>%
  group_by(WL_ORG) %>%
  mutate(
    hr_creatinine_imputed = if_else(WL_ORG == "HR", median_impute(hr_creatinine_for_score), as.numeric(hr_creatinine_for_score)),
    hr_bilirubin_imputed = if_else(WL_ORG == "HR", median_impute(hr_bilirubin), as.numeric(hr_bilirubin)),
    hr_albumin_imputed = if_else(WL_ORG == "HR", median_impute(hr_albumin), as.numeric(hr_albumin)),
    hr_sodium_imputed = if_else(WL_ORG == "HR", median_impute(hr_sodium), as.numeric(hr_sodium)),
    hr_bnp_imputed = if_else(WL_ORG == "HR", median_impute(hr_bnp), as.numeric(hr_bnp))
  ) %>%
  ungroup() %>%
  mutate(
    hr_egfr = calc_egfr_2021(hr_creatinine_imputed, age, sex),
    us_crs_proxy = -0.656 * hr_albumin_imputed +
      0.617 * safe_log1p(hr_bilirubin_imputed) -
      0.012 * hr_egfr -
      0.077 * hr_sodium_imputed -
      0.377 * hr_durable_lvad +
      1.092 * hr_short_mcs +
      0.433 * safe_log1p(hr_bnp_imputed),
    lung_las_cas_component_score = rowSums(cbind(
      as.numeric(coalesce(CAN_FEV1, 0)) * -0.01,
      as.numeric(coalesce(CAN_FVC, 0)) * -0.005,
      as.numeric(coalesce(CAN_PCO2, 0)) * 0.02,
      as.numeric(coalesce(CAN_AT_REST_O2, 0)) * 0.01,
      as.numeric(coalesce(CAN_SIX_MIN_WALK, 0)) * -0.002,
      as.numeric(coalesce(CAN_PULM_ART_MEAN, 0)) * 0.01,
      as.numeric(coalesce(CAN_CARDIAC_OUTPUT, 0)) * -0.05,
      flag_yes(coalesce(as.character(CAN_VENTILATOR), as.character(CAN_ON_VENTILATOR))) * 1.0,
      flag_yes(CAN_ECMO) * 1.5,
      flag_yes(CAN_CORTICOST_DEPND) * 0.3
    ), na.rm = TRUE),
    organ_score = case_when(
      WL_ORG == "HR" ~ us_crs_proxy,
      WL_ORG == "KI" ~ kidney_dialysis_age_score,
      WL_ORG == "LI" ~ liver_meld,
      WL_ORG == "LU" ~ lung_las_cas_component_score
    )
  ) %>%
  filter(
    index_year >= analysis_start_year,
    index_year <= analysis_end_year,
    CAN_SOURCE %in% c("A", "R"),
    !is.na(index_date),
    !is.na(observed_end_date),
    !is.na(followup_days),
    followup_days >= 0,
    observed_end_date >= index_date,
    !is.na(candidate_zip)
  ) %>%
  mutate(followup_days = pmax(followup_days, 0.5)) %>%
  select(
    waitlist_row_id, PERS_ID, PX_ID, WL_ORG, followup_days, adverse_event,
    age, sex, race, index_year_centered, organ_score, listing_center
  )

model_dat_all <- cohort %>%
  inner_join(exposures, by = c("waitlist_row_id", "PERS_ID", "PX_ID", "WL_ORG")) %>%
  filter(complete.cases(
    followup_days, adverse_event, age, sex, race, index_year_centered, organ_score, listing_center,
    tmean_waitlist_5c, temp_range_waitlist_5c, rhmean_waitlist_10pct,
    pm25_waitlist_5ug, o3_waitlist_10ppb, no2_waitlist_10unit
  ))

curve_specs <- tribble(
  ~exposure, ~scaled_term, ~raw_term, ~label, ~x_label,
  "Mean temperature", "tmean_waitlist_5c", "tmean_c_waitlist_mean", "Mean temperature", "Waitlist-period mean temperature (C)",
  "Diurnal temperature range", "temp_range_waitlist_5c", "temp_range_c_waitlist_mean", "Diurnal temperature range", "Waitlist-period mean diurnal temperature range (C)",
  "Mean relative humidity", "rhmean_waitlist_10pct", "rhmean_pct_waitlist_mean", "Mean relative humidity", "Waitlist-period mean relative humidity (%)"
)

base_adjusters <- c(
  "age", "sex", "race", "index_year_centered", "organ_score",
  "pm25_waitlist_5ug", "o3_waitlist_10ppb", "no2_waitlist_10unit",
  "strata(listing_center)", "cluster(PERS_ID)"
)

all_environment_terms <- c("tmean_waitlist_5c", "temp_range_waitlist_5c", "rhmean_waitlist_10pct")

make_curve <- function(dat, org, spec) {
  focal <- spec$scaled_term
  other_env <- setdiff(all_environment_terms, focal)
  rhs <- c(paste0("ns(", focal, ", df = 3)"), other_env, base_adjusters)
  form <- as.formula(paste("Surv(followup_days, adverse_event) ~", paste(rhs, collapse = " + ")))
  log_msg("Fitting ", org, " spline for ", spec$exposure, " n=", nrow(dat), " adverse=", sum(dat$adverse_event))
  fit <- coxph(form, data = dat, ties = "efron", robust = TRUE, x = FALSE, y = FALSE)

  x <- dat[[focal]]
  raw_x <- dat[[spec$raw_term]]
  qs <- quantile(x, probs = c(0.01, 0.50, 0.99), na.rm = TRUE)
  if (!all(is.finite(qs)) || qs[[1]] == qs[[3]]) return(tibble())

  grid <- seq(qs[[1]], qs[[3]], length.out = 120)
  raw_grid <- approx(x = x[order(x)], y = raw_x[order(x)], xout = grid, rule = 2, ties = mean)$y
  ref_value <- qs[[2]]
  ref_raw <- approx(x = x[order(x)], y = raw_x[order(x)], xout = ref_value, rule = 2, ties = mean)$y

  ref_idx <- which.min(abs(x - ref_value))
  ref_row <- dat[ref_idx, , drop = FALSE]
  newdat <- ref_row[rep(1, length(grid)), , drop = FALSE]
  newdat[[focal]] <- grid
  refdat <- ref_row
  refdat[[focal]] <- ref_value

  pred <- predict(fit, newdata = bind_rows(refdat, newdat), type = "lp", se.fit = TRUE)
  ref_fit <- pred$fit[[1]]
  curve_fit <- pred$fit[-1]
  curve_se <- pred$se.fit[-1]

  tibble(
    organ = org,
    organ_label = recode(org, !!!organ_labels),
    exposure = spec$exposure,
    label = spec$label,
    x_label = spec$x_label,
    exposure_value = raw_grid,
    reference_value = as.numeric(ref_raw),
    hazard_ratio = exp(curve_fit - ref_fit),
    conf_low = exp(curve_fit - ref_fit - 1.96 * curve_se),
    conf_high = exp(curve_fit - ref_fit + 1.96 * curve_se),
    n = nrow(dat),
    adverse_events = sum(dat$adverse_event),
    concordance = unname(summary(fit)$concordance[1])
  )
}

log_msg("Fitting exposure-response curves")
curves <- bind_rows(lapply(target_organs, function(org) {
  dat <- model_dat_all %>% filter(WL_ORG == org) %>% droplevels()
  bind_rows(lapply(seq_len(nrow(curve_specs)), function(i) make_curve(dat, org, curve_specs[i, ])))
}))

write_csv(curves, file.path(exposure_dir, "temperature_humidity_spline_response_curves.csv"))

theme_response <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )
}

plot_panel <- function(plot_dat, title, path) {
  p <- ggplot(plot_dat, aes(x = exposure_value, y = hazard_ratio)) +
    geom_hline(yintercept = 1, linewidth = 0.35, color = "grey45") +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.18, fill = "#3B6EA8") +
    geom_line(linewidth = 0.8, color = "#1E4F7A") +
    facet_wrap(~organ_label, ncol = 2, scales = "free_y") +
    scale_y_log10() +
    labs(
      title = title,
      subtitle = "Adjusted Cox spline curves; reference is organ-specific median exposure",
      x = unique(plot_dat$x_label),
      y = "Hazard ratio for death/deterioration delisting"
    ) +
    theme_response()

  ggsave(path, p, width = 11, height = 7, dpi = 240, bg = "white")
  path
}

paths <- bind_rows(lapply(unique(curves$exposure), function(exposure_name) {
  plot_dat <- curves %>% filter(exposure == exposure_name)
  path <- file.path(out_dir, paste0("exposure_response_", str_to_lower(str_replace_all(exposure_name, "[^A-Za-z0-9]+", "_")), ".png"))
  tibble(exposure = exposure_name, figure_path = plot_panel(plot_dat, paste0(exposure_name, " Exposure-Response by Organ"), path))
}))

combined_path <- file.path(out_dir, "temperature_humidity_exposure_response_all.png")
combined_curves <- curves %>%
  mutate(panel_label = paste(exposure, organ_label, sep = "\n"))
combined_plot <- ggplot(combined_curves, aes(x = exposure_value, y = hazard_ratio)) +
  geom_hline(yintercept = 1, linewidth = 0.3, color = "grey45") +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.16, fill = "#3B6EA8") +
  geom_line(linewidth = 0.7, color = "#1E4F7A") +
  facet_wrap(~panel_label, ncol = 4, scales = "free") +
  scale_y_log10() +
  labs(
    title = "Temperature and Humidity Exposure-Response by Organ",
    subtitle = "Adjusted Cox spline curves; reference is organ-specific median exposure",
    x = NULL,
    y = "Hazard ratio for death/deterioration delisting"
  ) +
  theme_response(base_size = 10) +
  theme(legend.position = "none")
ggsave(combined_path, combined_plot, width = 14, height = 9, dpi = 240, bg = "white")

paths <- bind_rows(paths, tibble(exposure = "All temperature/humidity", figure_path = combined_path))
write_csv(paths, file.path(out_dir, "temperature_humidity_response_figure_manifest.csv"))

log_msg("Wrote temperature/humidity response curves to ", out_dir)
