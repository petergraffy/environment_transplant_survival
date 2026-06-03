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

read_ozone_annual <- function() {
  read_csv(
    file.path(pollution_dir, "all_o3_zip.csv.gz"),
    col_types = cols(zip = col_character(), .default = col_guess()),
    show_col_types = FALSE
  ) %>%
    group_by(zip = clean_zip(zip), year) %>%
    summarise(o3_n_months = sum(!is.na(o3_ppb)), o3_ppb = mean(o3_ppb, na.rm = TRUE), .groups = "drop") %>%
    filter(o3_n_months >= 12L, !is.na(o3_ppb)) %>%
    mutate(o3_10ppb = o3_ppb / 10)
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
    lower_w <- if (lower == 0L) 0 else csum[[lower]]
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

log_msg("Reading annual ozone")
ozone_annual <- read_ozone_annual() %>%
  filter(year >= analysis_start_year - 1L, year <= analysis_end_year)

baseline_ozone <- ozone_annual %>%
  transmute(candidate_zip = zip, prior_exposure_year = year, o3_prior_1y_10ppb = o3_10ppb)

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
  left_join(baseline_ozone, by = c("candidate_zip", "prior_exposure_year")) %>%
  filter(complete.cases(kidney_dialysis_age_score, o3_prior_1y_10ppb))

log_msg("Constructing annual dynamic landmark rows")
intervals <- as.data.table(cohort)
intervals[, interval_year := Map(seq, index_year, as.integer(format(event_date, "%Y")))]
intervals <- intervals[, .(
  interval_year = unlist(interval_year)
), by = .(
  waitlist_row_id, PERS_ID, PX_ID, candidate_zip, index_date, end_date = event_date,
  followup_days, adverse_waitlist_event, kidney_dialysis_age_score, o3_prior_1y_10ppb
)]
intervals[, interval_start_date := pmax(index_date, as.Date(paste0(interval_year, "-01-01")))]
intervals[, interval_end_date := pmin(end_date, as.Date(paste0(interval_year, "-12-31")))]
intervals <- intervals[interval_end_date >= interval_start_date]
intervals[, landmark_days := as.numeric(interval_start_date - index_date)]
intervals[, tstop := pmax(as.numeric(end_date - interval_start_date), 0.5)]
intervals[, landmark_year := interval_year]
intervals[, exposure_year := landmark_year - 1L]

tv_exp <- as.data.table(ozone_annual %>%
  transmute(candidate_zip = zip, exposure_year = year, o3_current_10ppb = o3_10ppb))
setkey(intervals, candidate_zip, exposure_year)
setkey(tv_exp, candidate_zip, exposure_year)
intervals <- tv_exp[intervals]
intervals <- intervals[!is.na(o3_current_10ppb)]
setorder(intervals, waitlist_row_id, landmark_days)
intervals[, interval_days := pmax(as.numeric(interval_end_date - interval_start_date), 0.5)]
intervals[, o3_cumulative_10ppb := shift(cumsum(o3_current_10ppb * interval_days), fill = NA_real_) / shift(cumsum(interval_days), fill = NA_real_), by = waitlist_row_id]
intervals[is.na(o3_cumulative_10ppb), o3_cumulative_10ppb := o3_prior_1y_10ppb]

landmarks <- as_tibble(intervals) %>%
  mutate(
    landmark_year_factor = factor(landmark_year),
    landmark_days_scaled = landmark_days / 365.25
  ) %>%
  filter(landmark_days <= 5 * 365.25) %>%
  select(
    waitlist_row_id, PERS_ID, PX_ID, landmark_days, landmark_year, landmark_year_factor,
    landmark_days_scaled, tstop, adverse_waitlist_event,
    kidney_dialysis_age_score, o3_prior_1y_10ppb, o3_current_10ppb, o3_cumulative_10ppb
  ) %>%
  filter(complete.cases(kidney_dialysis_age_score, o3_prior_1y_10ppb, o3_current_10ppb, o3_cumulative_10ppb))

write_csv(
  landmarks %>%
    summarise(
      n_landmarks = n(),
      people = n_distinct(PERS_ID),
      adverse_events = sum(adverse_waitlist_event),
      median_landmark_days = median(landmark_days),
      p75_landmark_days = quantile(landmark_days, 0.75),
      cumulative_ozone_p25 = quantile(o3_cumulative_10ppb, 0.25),
      cumulative_ozone_p50 = quantile(o3_cumulative_10ppb, 0.50),
      cumulative_ozone_p75 = quantile(o3_cumulative_10ppb, 0.75)
    ),
  file.path(out_dir, "kidney_dynamic_landmark_ozone_auc_cohort_summary.csv")
)

set.seed(20260529)
people <- unique(landmarks$PERS_ID)
test_people <- sample(people, size = floor(length(people) * 0.30))
landmarks <- landmarks %>%
  mutate(split = if_else(PERS_ID %in% test_people, "test", "train"))
train <- landmarks %>% filter(split == "train")
test <- landmarks %>% filter(split == "test")

formulas <- list(
  dynamic_dialysis_age = Surv(tstop, adverse_waitlist_event) ~ kidney_dialysis_age_score + ns(landmark_days_scaled, df = 3),
  dynamic_dialysis_age_current_ozone = Surv(tstop, adverse_waitlist_event) ~ kidney_dialysis_age_score + ns(landmark_days_scaled, df = 3) + o3_current_10ppb,
  dynamic_dialysis_age_cumulative_ozone = Surv(tstop, adverse_waitlist_event) ~ kidney_dialysis_age_score + ns(landmark_days_scaled, df = 3) + o3_cumulative_10ppb,
  dynamic_dialysis_age_current_cumulative_ozone = Surv(tstop, adverse_waitlist_event) ~ kidney_dialysis_age_score + ns(landmark_days_scaled, df = 3) + o3_current_10ppb + o3_cumulative_10ppb
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
    model,
    train_n,
    test_n,
    train_events,
    test_events,
    coef = lapply(fit, function(x) {
      broom::tidy(x, exponentiate = TRUE, conf.int = TRUE) %>%
        select(term, estimate, conf.low, conf.high, p.value)
    })
  ) %>%
  tidyr::unnest(coef)

write_csv(auc_results, file.path(out_dir, "kidney_dynamic_landmark_ozone_time_dependent_auc.csv"))
write_csv(coef_results, file.path(out_dir, "kidney_dynamic_landmark_ozone_model_coefficients.csv"))
write_csv(predictions, file.path(out_dir, "kidney_dynamic_landmark_ozone_test_predictions.csv.gz"))

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
      dynamic_dialysis_age_current_ozone = "+ current ozone",
      dynamic_dialysis_age_cumulative_ozone = "+ cumulative ozone",
      dynamic_dialysis_age_current_cumulative_ozone = "+ current/cumulative ozone"
    ),
    horizon_label = factor(horizon_label, levels = unique(horizon_label))
  )

p_auc <- ggplot(plot_dat, aes(x = horizon_label, y = auc, color = model_label, group = model_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_y_continuous(labels = label_number(accuracy = 0.001), limits = c(0.50, NA)) +
  labs(
    title = "Kidney dynamic landmark prediction: test-set time-dependent AUC",
    x = "Prediction horizon after landmark",
    y = "IPCW time-dependent AUC",
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

p_delta <- ggplot(plot_dat %>% filter(model != "dynamic_dialysis_age"), aes(x = horizon_label, y = delta_auc_vs_dynamic_dialysis_age, color = model_label, group = model_label)) +
  geom_hline(yintercept = 0, color = "gray45", linewidth = 0.35) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_y_continuous(labels = label_number(accuracy = 0.0001)) +
  labs(
    title = "Incremental landmark AUC from updated ozone exposure",
    x = "Prediction horizon after landmark",
    y = "Delta AUC",
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "kidney_dynamic_landmark_ozone_time_dependent_auc.png"), p_auc, width = 8.5, height = 5.2, dpi = 320)
ggsave(file.path(fig_dir, "kidney_dynamic_landmark_ozone_incremental_auc.png"), p_delta, width = 8.5, height = 5.2, dpi = 320)

message("Wrote kidney dynamic landmark ozone AUC outputs to ", out_dir)
