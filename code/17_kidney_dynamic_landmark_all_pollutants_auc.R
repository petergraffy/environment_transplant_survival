#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(haven)
  library(readr)
  library(scales)
  library(splines)
  library(stringr)
  library(survival)
  library(tibble)
})

saf_paths_source <- file.path("code", "saf_paths.R")
if (!file.exists(saf_paths_source)) saf_paths_source <- "saf_paths.R"
source(saf_paths_source)
saf_paths <- get_saf_paths()
assert_saf_files(saf_paths, include_stathist = FALSE)

pubsaf_dir <- saf_paths$pubsaf_dir
supp_dir <- saf_paths$supp_dir
pollution_dir <- file.path("output", "zip_pollution")
out_dir <- file.path("output", "score_benchmark_pollution")
fig_dir <- file.path("output", "figures", "score_benchmark_pollution")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

analysis_start_year <- 2006L
analysis_end_year <- 2023L
horizons <- c(42, 90, 180, 365, 1095)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

clean_zip <- function(x) {
  z <- str_extract(as.character(x), "[0-9]+")
  z <- if_else(nchar(z) > 5L, substr(z, 1, 5), z)
  if_else(is.na(z), NA_character_, str_pad(z, 5, side = "left", pad = "0"))
}

flag_yes <- function(x) {
  y <- str_to_upper(str_trim(as.character(x)))
  as.integer(y %in% c("1", "Y", "YES", "TRUE", "T"))
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
    mutate(
      pm25_5ug = pm25_ug_m3 / 5,
      o3_10ppb = o3_ppb / 10,
      no2_10unit = no2 / 10
    )
}

ipcw_auc <- function(time, event, marker, horizon) {
  keep <- is.finite(time) & is.finite(event) & is.finite(marker) & time > 0
  time <- as.numeric(time[keep])
  event <- as.integer(event[keep])
  marker <- as.numeric(marker[keep])
  case_idx <- which(time <= horizon & event == 1L)
  control_idx <- which(time > horizon)
  if (length(case_idx) == 0L || length(control_idx) == 0L) {
    return(tibble(horizon_days = horizon, auc = NA_real_, n_cases = length(case_idx), n_controls = length(control_idx)))
  }

  censor_fit <- survfit(Surv(time, 1L - event) ~ 1)
  g_at <- function(t) {
    s <- summary(censor_fit, times = pmax(t, 0), extend = TRUE)$surv
    pmax(s, 1e-6)
  }
  case_w <- 1 / g_at(time[case_idx])
  control_w <- rep(1 / g_at(horizon), length(control_idx))
  control_marker <- marker[control_idx]
  ord <- order(control_marker)
  cm <- control_marker[ord]
  cw <- control_w[ord]
  csum <- cumsum(cw)
  total_control_w <- sum(cw)
  weighted_wins <- vapply(seq_along(case_idx), function(i) {
    m <- marker[case_idx[[i]]]
    lower <- findInterval(m, cm, left.open = FALSE)
    equal_left <- lower
    while (equal_left > 0L && cm[[equal_left]] == m) equal_left <- equal_left - 1L
    equal_right <- lower
    while (equal_right < length(cm) && cm[[equal_right + 1L]] == m) equal_right <- equal_right + 1L
    strictly_lower_w <- if (equal_left == 0L) 0 else csum[[equal_left]]
    equal_w <- if (equal_right >= equal_left + 1L) csum[[equal_right]] - strictly_lower_w else 0
    strictly_lower_w + 0.5 * equal_w
  }, numeric(1))
  auc <- sum(case_w * weighted_wins) / (sum(case_w) * total_control_w)
  tibble(horizon_days = horizon, auc = auc, n_cases = length(case_idx), n_controls = length(control_idx))
}

log_msg("Reading kidney ZIP supplement")
candidate_zip <- read_sas(file.path(supp_dir, "canzip2506.sas7bdat")) %>%
  transmute(PERS_ID, PX_ID, candidate_zip = clean_zip(CAN_PERM_ZIP)) %>%
  distinct(PERS_ID, PX_ID, .keep_all = TRUE)

log_msg("Reading annual pollution")
pollution_annual <- read_pollution_annual() %>%
  filter(year >= analysis_start_year - 1L, year <= analysis_end_year)

baseline_pollution <- pollution_annual %>%
  transmute(
    candidate_zip = zip,
    prior_exposure_year = year,
    pm25_prior_1y_5ug = pm25_5ug,
    o3_prior_1y_10ppb = o3_10ppb,
    no2_prior_1y_10unit = no2_10unit
  )

log_msg("Reading kidney candidate SAF")
kidney_cols <- c(
  "PERS_ID", "PX_ID", "WL_ORG", "CAN_LISTING_DT", "CAN_ACTIVATE_DT", "CAN_SOURCE",
  "CAN_REM_DT", "CAN_REM_CD", "CAN_DEATH_DT", "CAN_ENDWLFU", "CAN_AGE_AT_LISTING",
  "CAN_ON_DIAL", "CAN_DIAL", "CAN_DIAL_DT"
)

cohort <- read_sas(file.path(pubsaf_dir, "cand_kipa.sas7bdat"), col_select = any_of(kidney_cols)) %>%
  filter(WL_ORG == "KI") %>%
  left_join(candidate_zip, by = c("PERS_ID", "PX_ID")) %>%
  mutate(
    waitlist_row_id = row_number(),
    index_date = coalesce(CAN_ACTIVATE_DT, CAN_LISTING_DT),
    index_year = as.integer(format(index_date, "%Y")),
    prior_exposure_year = index_year - 1L,
    event_date = coalesce(CAN_REM_DT, CAN_DEATH_DT, CAN_ENDWLFU),
    followup_days = as.numeric(event_date - index_date),
    adverse_waitlist_event = as.integer(CAN_REM_CD %in% c(8, 13)),
    kidney_on_dialysis = flag_yes(coalesce(as.character(CAN_ON_DIAL), as.character(CAN_DIAL))),
    dialysis_years = as.numeric(index_date - CAN_DIAL_DT) / 365.25,
    dialysis_years = if_else(is.na(dialysis_years) | dialysis_years < 0, 0, dialysis_years),
    kidney_dialysis_age_score = CAN_AGE_AT_LISTING / 10 + 0.8 * kidney_on_dialysis + 0.3 * log1p(dialysis_years)
  ) %>%
  filter(
    index_year >= analysis_start_year,
    index_year <= analysis_end_year,
    CAN_SOURCE %in% c("A", "R"),
    !is.na(index_date),
    !is.na(event_date),
    !is.na(followup_days),
    followup_days >= 0,
    !is.na(candidate_zip)
  ) %>%
  mutate(followup_days = pmax(followup_days, 0.5)) %>%
  left_join(baseline_pollution, by = c("candidate_zip", "prior_exposure_year")) %>%
  filter(complete.cases(kidney_dialysis_age_score, pm25_prior_1y_5ug, o3_prior_1y_10ppb, no2_prior_1y_10unit))

log_msg("Constructing annual dynamic landmark rows")
intervals <- as.data.table(cohort)
intervals[, interval_year := Map(seq, index_year, as.integer(format(event_date, "%Y")))]
intervals <- intervals[, .(interval_year = unlist(interval_year)), by = .(
  waitlist_row_id, PERS_ID, PX_ID, candidate_zip, index_date, end_date = event_date,
  followup_days, adverse_waitlist_event, kidney_dialysis_age_score,
  pm25_prior_1y_5ug, o3_prior_1y_10ppb, no2_prior_1y_10unit
)]
intervals[, interval_start_date := pmax(index_date, as.Date(paste0(interval_year, "-01-01")))]
intervals[, interval_end_date := pmin(end_date, as.Date(paste0(interval_year, "-12-31")))]
intervals <- intervals[interval_end_date >= interval_start_date]
intervals[, landmark_days := as.numeric(interval_start_date - index_date)]
intervals[, tstop := pmax(as.numeric(end_date - interval_start_date), 0.5)]
intervals[, landmark_year := interval_year]
intervals[, exposure_year := landmark_year - 1L]

tv_exp <- as.data.table(pollution_annual %>%
  transmute(
    candidate_zip = zip,
    exposure_year = year,
    pm25_current_5ug = pm25_5ug,
    o3_current_10ppb = o3_10ppb,
    no2_current_10unit = no2_10unit
  ))
setkey(intervals, candidate_zip, exposure_year)
setkey(tv_exp, candidate_zip, exposure_year)
intervals <- tv_exp[intervals]
intervals <- intervals[complete.cases(pm25_current_5ug, o3_current_10ppb, no2_current_10unit)]
setorder(intervals, waitlist_row_id, landmark_days)
intervals[, interval_days := pmax(as.numeric(interval_end_date - interval_start_date), 0.5)]
intervals[, pm25_cumulative_5ug := shift(cumsum(pm25_current_5ug * interval_days), fill = NA_real_) / shift(cumsum(interval_days), fill = NA_real_), by = waitlist_row_id]
intervals[, o3_cumulative_10ppb := shift(cumsum(o3_current_10ppb * interval_days), fill = NA_real_) / shift(cumsum(interval_days), fill = NA_real_), by = waitlist_row_id]
intervals[, no2_cumulative_10unit := shift(cumsum(no2_current_10unit * interval_days), fill = NA_real_) / shift(cumsum(interval_days), fill = NA_real_), by = waitlist_row_id]
intervals[is.na(pm25_cumulative_5ug), pm25_cumulative_5ug := pm25_prior_1y_5ug]
intervals[is.na(o3_cumulative_10ppb), o3_cumulative_10ppb := o3_prior_1y_10ppb]
intervals[is.na(no2_cumulative_10unit), no2_cumulative_10unit := no2_prior_1y_10unit]

landmarks <- as_tibble(intervals) %>%
  mutate(landmark_days_scaled = landmark_days / 365.25) %>%
  filter(landmark_days <= 5 * 365.25) %>%
  select(
    waitlist_row_id, PERS_ID, PX_ID, landmark_days, landmark_year,
    landmark_days_scaled, tstop, adverse_waitlist_event,
    kidney_dialysis_age_score,
    pm25_current_5ug, o3_current_10ppb, no2_current_10unit,
    pm25_cumulative_5ug, o3_cumulative_10ppb, no2_cumulative_10unit
  ) %>%
  filter(complete.cases(.))

write_csv(
  landmarks %>%
    summarise(
      n_landmarks = n(),
      people = n_distinct(PERS_ID),
      adverse_events = sum(adverse_waitlist_event),
      median_landmark_days = median(landmark_days),
      p75_landmark_days = quantile(landmark_days, 0.75)
    ),
  file.path(out_dir, "kidney_dynamic_landmark_all_pollutants_auc_cohort_summary.csv")
)

set.seed(20260529)
people <- unique(landmarks$PERS_ID)
test_people <- sample(people, size = floor(length(people) * 0.30))
landmarks <- landmarks %>% mutate(split = if_else(PERS_ID %in% test_people, "test", "train"))
train <- landmarks %>% filter(split == "train")
test <- landmarks %>% filter(split == "test")

base_rhs <- "kidney_dialysis_age_score + ns(landmark_days_scaled, df = 3)"
formulas <- list(
  dynamic_dialysis_age = as.formula(paste("Surv(tstop, adverse_waitlist_event) ~", base_rhs)),
  current_pm25 = as.formula(paste("Surv(tstop, adverse_waitlist_event) ~", base_rhs, "+ pm25_current_5ug")),
  cumulative_pm25 = as.formula(paste("Surv(tstop, adverse_waitlist_event) ~", base_rhs, "+ pm25_cumulative_5ug")),
  current_o3 = as.formula(paste("Surv(tstop, adverse_waitlist_event) ~", base_rhs, "+ o3_current_10ppb")),
  cumulative_o3 = as.formula(paste("Surv(tstop, adverse_waitlist_event) ~", base_rhs, "+ o3_cumulative_10ppb")),
  current_no2 = as.formula(paste("Surv(tstop, adverse_waitlist_event) ~", base_rhs, "+ no2_current_10unit")),
  cumulative_no2 = as.formula(paste("Surv(tstop, adverse_waitlist_event) ~", base_rhs, "+ no2_cumulative_10unit")),
  current_all_pollutants = as.formula(paste("Surv(tstop, adverse_waitlist_event) ~", base_rhs, "+ pm25_current_5ug + o3_current_10ppb + no2_current_10unit")),
  cumulative_all_pollutants = as.formula(paste("Surv(tstop, adverse_waitlist_event) ~", base_rhs, "+ pm25_cumulative_5ug + o3_cumulative_10ppb + no2_cumulative_10unit"))
)

fits <- bind_rows(lapply(names(formulas), function(model_name) {
  log_msg("Fitting dynamic landmark model: ", model_name, " n_train=", nrow(train), " adverse=", sum(train$adverse_waitlist_event))
  fit <- coxph(formulas[[model_name]], data = train, ties = "efron", x = FALSE, y = FALSE)
  tibble(
    model = model_name,
    fit = list(fit),
    train_n = nrow(train),
    test_n = nrow(test),
    train_events = sum(train$adverse_waitlist_event),
    test_events = sum(test$adverse_waitlist_event)
  )
}))

predictions <- bind_rows(lapply(seq_len(nrow(fits)), function(i) {
  tibble(
    waitlist_row_id = test$waitlist_row_id,
    PERS_ID = test$PERS_ID,
    PX_ID = test$PX_ID,
    landmark_days = test$landmark_days,
    landmark_year = test$landmark_year,
    tstop = test$tstop,
    adverse_waitlist_event = test$adverse_waitlist_event,
    model = fits$model[[i]],
    marker = as.numeric(predict(fits$fit[[i]], newdata = test, type = "lp"))
  )
}))

auc_results <- predictions %>%
  group_by(model) %>%
  group_modify(~bind_rows(lapply(horizons, function(h) {
    ipcw_auc(.x$tstop, .x$adverse_waitlist_event, .x$marker, h)
  }))) %>%
  ungroup() %>%
  group_by(horizon_days) %>%
  mutate(
    baseline_auc = auc[model == "dynamic_dialysis_age"],
    delta_auc_vs_dynamic_dialysis_age = auc - baseline_auc
  ) %>%
  ungroup()

coef_results <- fits %>%
  transmute(
    model, train_n, test_n, train_events, test_events,
    coef = lapply(fit, function(x) {
      broom::tidy(x, exponentiate = TRUE, conf.int = TRUE) %>%
        select(term, estimate, conf.low, conf.high, p.value)
    })
  ) %>%
  tidyr::unnest(coef)

write_csv(auc_results, file.path(out_dir, "kidney_dynamic_landmark_all_pollutants_time_dependent_auc.csv"))
write_csv(coef_results, file.path(out_dir, "kidney_dynamic_landmark_all_pollutants_model_coefficients.csv"))
write_csv(predictions, file.path(out_dir, "kidney_dynamic_landmark_all_pollutants_test_predictions.csv.gz"))

plot_dat <- auc_results %>%
  mutate(
    horizon_label = case_when(
      horizon_days < 365 ~ paste0(horizon_days, " d"),
      horizon_days == 365 ~ "1 y",
      horizon_days == 1095 ~ "3 y",
      TRUE ~ paste0(round(horizon_days / 365.25, 1), " y")
    ),
    model_label = recode(
      model,
      dynamic_dialysis_age = "Dynamic dialysis + age",
      current_pm25 = "+ current PM2.5",
      cumulative_pm25 = "+ cumulative PM2.5",
      current_o3 = "+ current ozone",
      cumulative_o3 = "+ cumulative ozone",
      current_no2 = "+ current NO2",
      cumulative_no2 = "+ cumulative NO2",
      current_all_pollutants = "+ current all pollutants",
      cumulative_all_pollutants = "+ cumulative all pollutants"
    ),
    horizon_label = factor(horizon_label, levels = unique(horizon_label))
  )

p_delta <- ggplot(plot_dat %>% filter(model != "dynamic_dialysis_age"), aes(x = horizon_label, y = delta_auc_vs_dynamic_dialysis_age, color = model_label, group = model_label)) +
  geom_hline(yintercept = 0, color = "gray45", linewidth = 0.35) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.8) +
  scale_y_continuous(labels = label_number(accuracy = 0.0001)) +
  labs(
    title = "Incremental kidney dynamic landmark AUC from updated pollution",
    x = "Prediction horizon after landmark",
    y = "Delta AUC",
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "kidney_dynamic_landmark_all_pollutants_incremental_auc.png"), p_delta, width = 9.8, height = 6.0, dpi = 320)

message("Wrote kidney dynamic landmark all-pollutant AUC outputs to ", out_dir)
