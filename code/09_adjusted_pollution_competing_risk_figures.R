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
  library(nnet)
  library(readr)
  library(scales)
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
assert_saf_files(saf_paths)

pubsaf_dir <- saf_paths$pubsaf_dir
supp_dir <- saf_paths$supp_dir
pollution_dir <- file.path("output", "zip_pollution")
community_dir <- file.path("data", "processed", "community")
out_dir <- file.path("output", "figures", "adjusted_competing_risks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analysis_start_year <- 2006L
analysis_end_year <- 2023L
target_organs <- c("HR", "KI", "LI", "LU")
max_curve_days <- 5 * 365.25
curve_grid <- seq(0, max_curve_days, by = 30.4375)
max_model_n <- as.integer(Sys.getenv("ADJUSTED_CIF_MODEL_N", "100000"))
max_standardize_n <- as.integer(Sys.getenv("ADJUSTED_CIF_STANDARDIZE_N", "5000"))
max_ps_n <- as.integer(Sys.getenv("ADJUSTED_CIF_PS_N", "120000"))

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
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

scale_01_to_pct <- function(x) {
  as.numeric(x) * 100
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
  female <- str_to_upper(as.character(sex)) %in% c("F", "FEMALE")
  kappa <- if_else(female, 0.7, 0.9)
  alpha <- if_else(female, -0.241, -0.302)
  sex_mult <- if_else(female, 1.012, 1.0)
  142 * pmin(scr / kappa, 1)^alpha * pmax(scr / kappa, 1)^(-1.200) * 0.9938^age * sex_mult
}

read_pollution_annual <- function() {
  col_types <- cols(zip = col_character(), .default = col_guess())
  pm25 <- read_csv(file.path(pollution_dir, "all_pm25_zip.csv.gz"), col_types = col_types, show_col_types = FALSE) %>%
    group_by(zip = clean_zip(zip), year) %>%
    summarise(pm25_n_months = sum(!is.na(pm25_ug_m3)), pm25_ug_m3 = mean(pm25_ug_m3, na.rm = TRUE), .groups = "drop")
  o3 <- read_csv(file.path(pollution_dir, "all_o3_zip.csv.gz"), col_types = col_types, show_col_types = FALSE) %>%
    group_by(zip = clean_zip(zip), year) %>%
    summarise(o3_n_months = sum(!is.na(o3_ppb)), o3_ppb = mean(o3_ppb, na.rm = TRUE), .groups = "drop")
  no2 <- read_csv(file.path(pollution_dir, "all_no2_zip.csv.gz"), col_types = col_types, show_col_types = FALSE) %>%
    transmute(zip = clean_zip(zip), year, no2 = no2)

  pm25 %>%
    inner_join(o3, by = c("zip", "year")) %>%
    inner_join(no2, by = c("zip", "year")) %>%
    filter(pm25_n_months >= 12L, o3_n_months >= 12L, !is.na(pm25_ug_m3), !is.na(o3_ppb), !is.na(no2)) %>%
    mutate(pm25_5ug = pm25_ug_m3 / 5, o3_10ppb = o3_ppb / 10, no2_10unit = no2 / 10)
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
  "CAN_GENDER", "CAN_RACE_SRTR", "CAN_ETHNICITY_SRTR", "CAN_ABO", "CAN_LISTING_CTR_CD",
  "CAN_ON_DIAL", "CAN_DIAL", "CAN_DIAL_DT", "CAN_MOST_RECENT_CREAT", "CAN_TOT_ALBUMIN",
  "CAN_TOT_BILI", "CAN_LAST_SERUM_SODIUM", "CAN_INIT_SRTR_LAB_MELD", "CAN_LAST_SRTR_LAB_MELD",
  "CAN_FVC", "CAN_FEV1", "CAN_PCO2", "CAN_AT_REST_O2", "CAN_SIX_MIN_WALK",
  "CAN_PULM_ART_MEAN", "CAN_CARDIAC_OUTPUT", "CAN_VENTILATOR", "CAN_ON_VENTILATOR",
  "CAN_ECMO", "CAN_CORTICOST_DEPND"
)

log_msg("Reading inputs")
candidate_zip <- read_sas(file.path(supp_dir, "canzip2506.sas7bdat")) %>%
  transmute(PERS_ID, PX_ID, candidate_zip = clean_zip(CAN_PERM_ZIP)) %>%
  distinct(PERS_ID, PX_ID, .keep_all = TRUE)

pollution_annual <- read_pollution_annual() %>%
  filter(year >= analysis_start_year - 1L, year <= analysis_end_year)

community <- read_csv(
  file.path(community_dir, "zcta_acs_community_covariates_2005_2023.csv.gz"),
  col_types = cols(zip = col_character(), .default = col_guess()),
  show_col_types = FALSE
) %>%
  transmute(
    zip = clean_zip(zip),
    analysis_year,
    median_household_income,
    pct_poverty = scale_01_to_pct(pct_poverty),
    pct_bachelor_plus = scale_01_to_pct(pct_bachelor_plus),
    pct_unemployed = scale_01_to_pct(pct_unemployed),
    pct_no_vehicle = scale_01_to_pct(pct_no_vehicle),
    pct_black = scale_01_to_pct(pct_black),
    pct_hispanic = scale_01_to_pct(pct_hispanic)
  )

candidate_data <- candidate_files %>%
  rowwise() %>%
  mutate(data = list(read_sas(path, col_select = any_of(candidate_cols)) %>% mutate(candidate_group = candidate_group))) %>%
  ungroup() %>%
  select(data) %>%
  unnest(data)

heart_just <- read_sas(
  file.path(pubsaf_dir, "statjust_hr1a.sas7bdat"),
  col_select = any_of(c(
    "PX_ID", "WL_ORG", "CANHX_CHG_DT", "CANHX_VAD", "CANHX_LVAD_TYPE", "CANHX_RVAD_TYPE",
    "CANHX_ECMO", "CANHX_LAB_SERUM_CREAT", "CANHX_LAB_BILI", "CANHX_LAB_ALBUMIN",
    "CANHX_LAB_SODIUM", "CANHX_LAB_BNP"
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

baseline_exposure <- pollution_annual %>%
  transmute(
    candidate_zip = zip,
    prior_exposure_year = year,
    pm25_ug_m3 = pm25_ug_m3,
    o3_ppb = o3_ppb,
    no2 = no2
  )

log_msg("Constructing analysis cohort")
cohort <- candidate_data %>%
  filter(WL_ORG %in% target_organs) %>%
  left_join(candidate_zip, by = c("PERS_ID", "PX_ID")) %>%
  left_join(heart_just, by = "PX_ID") %>%
  mutate(
    index_date = coalesce(CAN_ACTIVATE_DT, CAN_LISTING_DT),
    index_year = as.integer(format(index_date, "%Y")),
    prior_exposure_year = index_year - 1L,
    event_date = coalesce(CAN_REM_DT, CAN_DEATH_DT, CAN_ENDWLFU),
    followup_days = as.numeric(event_date - index_date),
    adverse_event = as.integer(CAN_REM_CD %in% c(8, 13)),
    transplant_event = as.integer(CAN_REM_CD %in% c(4, 14, 15, 18, 19, 21, 22, 23)),
    event_type = case_when(adverse_event == 1L ~ 1L, transplant_event == 1L ~ 2L, TRUE ~ 0L),
    age = CAN_AGE_AT_LISTING,
    sex = as_factor_missing(if_else(CAN_GENDER %in% c("M", "F"), CAN_GENDER, NA_character_)),
    race_srtr = as_factor_missing(CAN_RACE_SRTR),
    ethnicity_srtr = as_factor_missing(CAN_ETHNICITY_SRTR),
    listing_center = collapse_rare(CAN_LISTING_CTR_CD, min_n = 100L),
    index_year_factor = factor(index_year),
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
    hr_durable_lvad = coalesce(hr_durable_lvad, 0L)
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
    benchmark_score = case_when(
      WL_ORG == "HR" ~ us_crs_proxy,
      WL_ORG == "KI" ~ kidney_dialysis_age_score,
      WL_ORG == "LI" ~ liver_meld,
      WL_ORG == "LU" ~ lung_las_cas_component_score
    ),
    organ = factor(WL_ORG, levels = target_organs, labels = c("Heart", "Kidney", "Liver", "Lung"))
  ) %>%
  filter(
    index_year >= analysis_start_year,
    index_year <= analysis_end_year,
    CAN_SOURCE %in% c("A", "R"),
    !is.na(index_date), !is.na(event_date), !is.na(followup_days),
    followup_days >= 0, !is.na(candidate_zip)
  ) %>%
  mutate(followup_days = pmax(followup_days, 0.5)) %>%
  left_join(baseline_exposure, by = c("candidate_zip", "prior_exposure_year")) %>%
  left_join(community, by = c("candidate_zip" = "zip", "index_year" = "analysis_year")) %>%
  filter(!is.na(benchmark_score), !is.na(pm25_ug_m3), !is.na(o3_ppb), !is.na(no2))

write_csv(
  cohort %>%
    group_by(organ) %>%
    summarise(n = n(), adverse_events = sum(event_type == 1L), transplants = sum(event_type == 2L), .groups = "drop"),
  file.path(out_dir, "adjusted_competing_risk_cohort_summary.csv")
)

community_terms <- c(
  "ns(median_household_income, df = 3)", "ns(pct_poverty, df = 3)",
  "ns(pct_bachelor_plus, df = 3)", "ns(pct_unemployed, df = 3)",
  "ns(pct_no_vehicle, df = 3)", "ns(pct_black, df = 3)", "ns(pct_hispanic, df = 3)"
)

weighted_cif <- function(dat, qlevel, max_time, grid) {
  qdat <- dat %>%
    filter(exposure_quartile == qlevel, followup_days <= max_time | followup_days > 0) %>%
    mutate(time = pmin(followup_days, max_time), event_at_time = if_else(followup_days <= max_time, event_type, 0L))

  if (nrow(qdat) == 0L) return(tibble())

  tab <- qdat %>%
    group_by(time) %>%
    summarise(
      removed_w = sum(iptw, na.rm = TRUE),
      adverse_w = sum(iptw[event_at_time == 1L], na.rm = TRUE),
      transplant_w = sum(iptw[event_at_time == 2L], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(time) %>%
    mutate(risk_w = rev(cumsum(rev(removed_w))))

  surv <- 1
  cif1 <- 0
  cif2 <- 0
  out <- vector("list", nrow(tab))
  for (i in seq_len(nrow(tab))) {
    if (tab$risk_w[[i]] > 0) {
      h1 <- tab$adverse_w[[i]] / tab$risk_w[[i]]
      h2 <- tab$transplant_w[[i]] / tab$risk_w[[i]]
      cif1 <- cif1 + surv * h1
      cif2 <- cif2 + surv * h2
      surv <- surv * pmax(0, 1 - h1 - h2)
    }
    out[[i]] <- tibble(time_days = tab$time[[i]], adverse = cif1, transplant = cif2)
  }

  bind_rows(out) %>%
    complete(time_days = unique(c(0, time_days, grid))) %>%
    arrange(time_days) %>%
    fill(adverse, transplant, .direction = "down") %>%
    mutate(adverse = replace_na(adverse, 0), transplant = replace_na(transplant, 0)) %>%
    filter(time_days %in% grid) %>%
    pivot_longer(c(adverse, transplant), names_to = "event", values_to = "cumulative_incidence") %>%
    mutate(
      event = recode(event, adverse = "Death/delisted too sick", transplant = "Transplant"),
      quartile = qlevel
    )
}

make_pollutant_curves <- function(pollutant, value_col, label, unit) {
  log_msg("Adjusted curves for ", label)
  by_organ <- cohort %>%
    group_by(organ) %>%
    mutate(
      exposure_quartile = ntile(.data[[value_col]], 4),
      exposure_quartile = factor(paste0("Q", exposure_quartile), levels = paste0("Q", 1:4))
    ) %>%
    ungroup()

  quartile_labels <- by_organ %>%
    group_by(organ, exposure_quartile) %>%
    summarise(
      min_value = min(.data[[value_col]], na.rm = TRUE),
      max_value = max(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      pollutant = label,
      quartile_label = paste0(exposure_quartile, " (", round(min_value, 1), "-", round(max_value, 1), " ", unit, ")")
    )

  curves <- lapply(levels(by_organ$organ), function(org) {
    dat <- by_organ %>% filter(organ == org)
    if (nrow(dat) == 0 || sum(dat$event_type == 1L) < 100 || sum(dat$event_type == 2L) < 100) return(NULL)

    model_dat <- dat %>%
      filter(complete.cases(
        followup_days, event_type, exposure_quartile, benchmark_score, age, sex, race_srtr,
        ethnicity_srtr, index_year_factor, median_household_income, pct_poverty,
        pct_bachelor_plus, pct_unemployed, pct_no_vehicle, pct_black, pct_hispanic
      )) %>%
      droplevels()

    ps_formula <- exposure_quartile ~ ns(benchmark_score, df = 3) + ns(age, df = 4) +
      sex + race_srtr + ethnicity_srtr + index_year_factor +
      ns(median_household_income, df = 3) + ns(pct_poverty, df = 3) +
      ns(pct_bachelor_plus, df = 3) + ns(pct_unemployed, df = 3) +
      ns(pct_no_vehicle, df = 3) + ns(pct_black, df = 3) + ns(pct_hispanic, df = 3)

    if (nrow(model_dat) > max_ps_n) {
      set.seed(20260528 + 100 * match(label, c("PM2.5", "Ozone", "NO2")) + match(org, levels(by_organ$organ)))
      ps_fit_dat <- model_dat %>%
        group_by(exposure_quartile, event_type) %>%
        slice_sample(prop = min(1, max_ps_n / nrow(model_dat))) %>%
        ungroup() %>%
        droplevels()
    } else {
      ps_fit_dat <- model_dat
    }

    log_msg("Fitting IPTW exposure model for ", label, " ", org, " ps_n=", nrow(ps_fit_dat), " full_n=", nrow(model_dat))
    ps_fit <- multinom(ps_formula, data = ps_fit_dat, trace = FALSE, MaxNWts = 20000, maxit = 200)
    ps <- predict(ps_fit, newdata = model_dat, type = "probs")
    if (is.vector(ps)) ps <- cbind(Q1 = 1 - ps, Q2 = ps)
    ps <- as.matrix(ps[, levels(model_dat$exposure_quartile), drop = FALSE])
    marginal <- prop.table(table(model_dat$exposure_quartile))
    q_index <- cbind(seq_len(nrow(model_dat)), match(model_dat$exposure_quartile, colnames(ps)))
    model_dat$iptw <- as.numeric(marginal[as.character(model_dat$exposure_quartile)]) / pmax(ps[q_index], 0.01)
    trim <- quantile(model_dat$iptw, probs = c(0.01, 0.99), na.rm = TRUE)
    model_dat$iptw <- pmin(pmax(model_dat$iptw, trim[[1]]), trim[[2]])

    bind_rows(lapply(levels(model_dat$exposure_quartile), function(q) weighted_cif(model_dat, q, max_curve_days, curve_grid))) %>%
      mutate(organ = org, pollutant = label, ps_n = nrow(ps_fit_dat), full_n = nrow(model_dat))
  }) %>%
    bind_rows() %>%
    left_join(quartile_labels, by = c("organ", "quartile" = "exposure_quartile", "pollutant")) %>%
    mutate(
      years_since_listing = time_days / 365.25,
      quartile = factor(quartile, levels = paste0("Q", 1:4))
    )

  write_csv(curves, file.path(out_dir, paste0("adjusted_competing_risk_curves_", pollutant, ".csv")))
  write_csv(quartile_labels, file.path(out_dir, paste0("adjusted_competing_risk_quartiles_", pollutant, ".csv")))

  p <- ggplot(curves, aes(x = years_since_listing, y = cumulative_incidence, color = quartile)) +
    geom_line(linewidth = 0.8) +
    facet_grid(organ ~ event, scales = "free_y") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_color_brewer(palette = "Dark2") +
    labs(
      title = paste0("Adjusted competing risks by ", label, " quartile"),
      subtitle = "IPTW-adjusted Aalen-Johansen CIFs; adverse outcome is death or delisted too sick",
      x = "Years since listing",
      y = "Cumulative incidence",
      color = "Quartile"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggsave(file.path(out_dir, paste0("adjusted_competing_risks_", pollutant, "_quartile.png")), p, width = 11, height = 9, dpi = 300)
  curves
}

all_curves <- bind_rows(
  make_pollutant_curves("pm25", "pm25_ug_m3", "PM2.5", "ug/m3"),
  make_pollutant_curves("o3", "o3_ppb", "Ozone", "ppb"),
  make_pollutant_curves("no2", "no2", "NO2", "units")
)

write_csv(all_curves, file.path(out_dir, "adjusted_competing_risk_curves_all_pollutants.csv"))
message("Wrote adjusted competing-risk figures to ", out_dir)
