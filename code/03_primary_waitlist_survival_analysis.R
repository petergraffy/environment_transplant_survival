#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(ggplot2)
  library(haven)
  library(readr)
  library(splines)
  library(stringr)
  library(survival)
  library(tibble)
})

saf_paths_source <- file.path("code", "saf_paths.R")
if (!file.exists(saf_paths_source)) saf_paths_source <- "saf_paths.R"
source(saf_paths_source)
saf_paths <- get_saf_paths()
assert_saf_files(saf_paths)
saf_dir <- saf_paths$saf_dir
pubsaf_dir <- saf_paths$pubsaf_dir
supp_dir <- saf_paths$supp_dir
pollution_dir <- "output/zip_pollution"
out_dir <- "output/primary_survival"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Using SAF directory: ", saf_dir)

analysis_start_year <- 2005L
analysis_end_year <- 2023L

candidate_files <- tribble(
  ~candidate_group, ~path,
  "kidney_pancreas", file.path(pubsaf_dir, "cand_kipa.sas7bdat"),
  "liver_intestine", file.path(pubsaf_dir, "cand_liin.sas7bdat"),
  "thoracic", file.path(pubsaf_dir, "cand_thor.sas7bdat")
)

clean_zip <- function(x) {
  z <- str_extract(as.character(x), "\\d+")
  z <- if_else(nchar(z) > 5L, substr(z, 1, 5), z)
  if_else(is.na(z), NA_character_, str_pad(z, 5, side = "left", pad = "0"))
}

mode_value <- function(x) {
  ux <- na.omit(x)
  if (length(ux) == 0L) return(NA_character_)
  names(which.max(table(ux)))
}

message("Reading candidate ZIP supplement")
candidate_zip <- read_sas(file.path(supp_dir, "canzip2506.sas7bdat")) %>%
  transmute(
    PERS_ID,
    PX_ID,
    candidate_zip = clean_zip(CAN_PERM_ZIP)
  ) %>%
  distinct(PERS_ID, PX_ID, .keep_all = TRUE)

message("Building annual pollution exposure table")
pollution_col_types <- cols(zip = col_character(), .default = col_guess())

pm25_year <- read_csv(
  file.path(pollution_dir, "all_pm25_zip.csv.gz"),
  show_col_types = FALSE,
  col_types = pollution_col_types
) %>%
  filter(year >= analysis_start_year, year <= analysis_end_year) %>%
  group_by(zip = clean_zip(zip), year) %>%
  summarise(
    pm25_n_months = sum(!is.na(pm25_ug_m3)),
    pm25_ug_m3 = mean(pm25_ug_m3, na.rm = TRUE),
    pm25_source = mode_value(value_source),
    .groups = "drop"
  )

o3_year <- read_csv(
  file.path(pollution_dir, "all_o3_zip.csv.gz"),
  show_col_types = FALSE,
  col_types = pollution_col_types
) %>%
  filter(year >= analysis_start_year, year <= analysis_end_year) %>%
  group_by(zip = clean_zip(zip), year) %>%
  summarise(
    o3_n_months = sum(!is.na(o3_ppb)),
    o3_ppb = mean(o3_ppb, na.rm = TRUE),
    o3_source = mode_value(value_source),
    .groups = "drop"
  )

no2_year <- read_csv(
  file.path(pollution_dir, "all_no2_zip.csv.gz"),
  show_col_types = FALSE,
  col_types = pollution_col_types
) %>%
  filter(year >= analysis_start_year, year <= analysis_end_year) %>%
  transmute(
    zip = clean_zip(zip),
    year,
    no2 = no2,
    no2_source = value_source
  )

write_csv(
  bind_rows(
    pm25_year %>% summarise(table = "pm25_year", n_rows = n(), n_zips = n_distinct(zip), min_year = min(year), max_year = max(year)),
    o3_year %>% summarise(table = "o3_year", n_rows = n(), n_zips = n_distinct(zip), min_year = min(year), max_year = max(year)),
    no2_year %>% summarise(table = "no2_year", n_rows = n(), n_zips = n_distinct(zip), min_year = min(year), max_year = max(year)),
    inner_join(pm25_year, o3_year, by = c("zip", "year")) %>%
      summarise(table = "pm25_o3_join", n_rows = n(), n_zips = n_distinct(zip), min_year = min(year), max_year = max(year)),
    inner_join(pm25_year, no2_year, by = c("zip", "year")) %>%
      summarise(table = "pm25_no2_join", n_rows = n(), n_zips = n_distinct(zip), min_year = min(year), max_year = max(year))
  ),
  file.path(out_dir, "annual_pollution_exposure_debug.csv")
)

pollution_year <- pm25_year %>%
  inner_join(o3_year, by = c("zip", "year")) %>%
  inner_join(no2_year, by = c("zip", "year")) %>%
  filter(
    !is.na(pm25_ug_m3),
    !is.na(o3_ppb),
    !is.na(no2),
    pm25_n_months >= 12L,
    o3_n_months >= 12L
  ) %>%
  mutate(
    pm25_5ug = pm25_ug_m3 / 5,
    o3_10ppb = o3_ppb / 10,
    no2_10unit = no2 / 10,
    any_nearest_fill = pm25_source != "zonal_mean" | o3_source != "zonal_mean" | no2_source != "zonal_mean"
  )

if (nrow(pollution_year) == 0L) {
  stop("Annual pollution exposure table is empty; see annual_pollution_exposure_debug.csv")
}

write_csv(
  pollution_year %>%
    summarise(
      n_zip_years = n(),
      n_zips = n_distinct(zip),
      min_year = min(year),
      max_year = max(year)
    ),
  file.path(out_dir, "annual_pollution_exposure_summary.csv")
)

message("Reading candidate SAF files")
candidate_data <- candidate_files %>%
  rowwise() %>%
  mutate(data = list(
    read_sas(
      path,
      col_select = any_of(c(
        "PERS_ID", "PX_ID", "WL_ORG", "CAN_LISTING_DT", "CAN_ACTIVATE_DT",
        "CAN_SOURCE", "CAN_REM_DT", "CAN_REM_CD", "CAN_DEATH_DT", "CAN_ENDWLFU",
        "CAN_AGE_AT_LISTING", "CAN_GENDER", "CAN_RACE_SRTR", "CAN_ETHNICITY_SRTR",
        "CAN_BMI", "CAN_MED_COND", "CAN_LISTING_CTR_CD"
      ))
    ) %>%
      mutate(candidate_group = candidate_group)
  )) %>%
  ungroup() %>%
  select(data) %>%
  tidyr::unnest(data)

message("Constructing survival cohort")
cohort <- candidate_data %>%
  left_join(candidate_zip, by = c("PERS_ID", "PX_ID")) %>%
  mutate(
    index_date = coalesce(CAN_ACTIVATE_DT, CAN_LISTING_DT),
    index_year = as.integer(format(index_date, "%Y")),
    event_waitlist_death = CAN_REM_CD == 8,
    event_date = if_else(event_waitlist_death, coalesce(CAN_DEATH_DT, CAN_REM_DT), as.Date(NA)),
    censor_date = coalesce(CAN_REM_DT, CAN_ENDWLFU),
    end_date = if_else(event_waitlist_death, event_date, censor_date),
    followup_days = as.numeric(end_date - index_date),
    age = CAN_AGE_AT_LISTING,
    sex = factor(if_else(CAN_GENDER %in% c("M", "F"), CAN_GENDER, NA_character_)),
    race_srtr = factor(na_if(CAN_RACE_SRTR, "")),
    ethnicity_srtr = factor(na_if(CAN_ETHNICITY_SRTR, "")),
    wl_org = factor(WL_ORG),
    medical_condition = factor(CAN_MED_COND),
    listing_center = factor(na_if(CAN_LISTING_CTR_CD, ""))
  ) %>%
  filter(
    index_year >= analysis_start_year,
    index_year <= analysis_end_year,
    CAN_SOURCE %in% c("A", "R"),
    !is.na(index_date),
    !is.na(end_date),
    !is.na(followup_days),
    followup_days >= 0
  ) %>%
  mutate(
    followup_days = pmax(followup_days, 0.5),
    event_waitlist_death = as.integer(event_waitlist_death)
  ) %>%
  left_join(
    pollution_year,
    by = c("candidate_zip" = "zip", "index_year" = "year")
  ) %>%
  mutate(
    complete_pollution = !is.na(pm25_ug_m3) & !is.na(o3_ppb) & !is.na(no2),
    complete_model_covariates = complete.cases(
      age, sex, race_srtr, ethnicity_srtr, CAN_BMI, wl_org, index_year
    )
  )

analysis_cohort <- cohort %>%
  filter(complete_pollution, complete_model_covariates) %>%
  mutate(index_year_factor = factor(index_year))

write_csv(
  cohort %>%
    summarise(
      n_registrations = n(),
      n_people = n_distinct(PERS_ID),
      n_waitlist_deaths = sum(event_waitlist_death),
      n_complete_pollution = sum(complete_pollution),
      n_complete_covariates = sum(complete_model_covariates),
      n_primary_analysis = sum(complete_pollution & complete_model_covariates)
    ),
  file.path(out_dir, "primary_survival_cohort_summary.csv")
)

write_csv(
  analysis_cohort %>%
    group_by(WL_ORG) %>%
    summarise(
      n_registrations = n(),
      n_people = n_distinct(PERS_ID),
      n_waitlist_deaths = sum(event_waitlist_death),
      median_followup_days = median(followup_days),
      mean_pm25_ug_m3 = mean(pm25_ug_m3),
      mean_o3_ppb = mean(o3_ppb),
      mean_no2 = mean(no2),
      .groups = "drop"
    ) %>%
    arrange(WL_ORG),
  file.path(out_dir, "primary_survival_by_organ_summary.csv")
)

write_csv(
  analysis_cohort %>%
    summarise(
      n = n(),
      cor_pm25_o3 = cor(pm25_ug_m3, o3_ppb),
      cor_pm25_no2 = cor(pm25_ug_m3, no2),
      cor_o3_no2 = cor(o3_ppb, no2)
    ),
  file.path(out_dir, "pollutant_correlation_summary.csv")
)

fit_cox <- function(formula, data, model_name) {
  message("Fitting ", model_name)
  fit <- coxph(
    formula,
    data = data,
    ties = "efron",
    robust = TRUE,
    x = FALSE,
    y = FALSE
  )
  list(
    fit = fit,
    table = tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(model = model_name, .before = 1)
  )
}

base_adjustment <- Surv(followup_days, event_waitlist_death) ~
  ns(age, df = 4) + sex + race_srtr + ethnicity_srtr + CAN_BMI +
  strata(wl_org) + strata(index_year_factor) + cluster(PERS_ID)

models <- list(
  fit_cox(
    update(base_adjustment, . ~ . + pm25_5ug),
    analysis_cohort,
    "single_pollutant_pm25_per_5ug_m3"
  ),
  fit_cox(
    update(base_adjustment, . ~ . + o3_10ppb),
    analysis_cohort,
    "single_pollutant_o3_per_10ppb"
  ),
  fit_cox(
    update(base_adjustment, . ~ . + no2_10unit),
    analysis_cohort,
    "single_pollutant_no2_per_10unit"
  ),
  fit_cox(
    update(base_adjustment, . ~ . + pm25_5ug + o3_10ppb + no2_10unit),
    analysis_cohort,
    "three_pollutant_model"
  )
)

cox_results <- bind_rows(lapply(models, `[[`, "table"))
write_csv(cox_results, file.path(out_dir, "cox_model_results_full.csv"))

pollutant_results <- cox_results %>%
  filter(term %in% c("pm25_5ug", "o3_10ppb", "no2_10unit")) %>%
  transmute(
    model,
    pollutant = recode(
      term,
      pm25_5ug = "PM2.5 per 5 ug/m3",
      o3_10ppb = "O3 per 10 ppb",
      no2_10unit = "NO2 per 10 units"
    ),
    hazard_ratio = estimate,
    conf_low = conf.low,
    conf_high = conf.high,
    robust_se = std.error,
    statistic,
    p_value = p.value
  )
write_csv(pollutant_results, file.path(out_dir, "cox_pollutant_results.csv"))

organ_results <- analysis_cohort %>%
  group_split(WL_ORG) %>%
  lapply(function(dat) {
    org <- unique(dat$WL_ORG)
    if (sum(dat$event_waitlist_death) < 100L) return(NULL)
    form <- Surv(followup_days, event_waitlist_death) ~
      pm25_5ug + o3_10ppb + no2_10unit +
      ns(age, df = 4) + sex + race_srtr + ethnicity_srtr + CAN_BMI +
      strata(index_year_factor) + cluster(PERS_ID)
    fit <- coxph(form, data = dat, ties = "efron", robust = TRUE, x = FALSE, y = FALSE)
    tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term %in% c("pm25_5ug", "o3_10ppb", "no2_10unit")) %>%
      mutate(WL_ORG = as.character(org), n = nrow(dat), deaths = sum(dat$event_waitlist_death), .before = 1)
  }) %>%
  bind_rows()

write_csv(organ_results, file.path(out_dir, "cox_pollutant_results_by_organ.csv"))

plot_data <- pollutant_results %>%
  filter(model == "three_pollutant_model") %>%
  mutate(pollutant = factor(pollutant, levels = rev(c(
    "PM2.5 per 5 ug/m3", "O3 per 10 ppb", "NO2 per 10 units"
  ))))

ggplot(plot_data, aes(x = hazard_ratio, y = pollutant)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.15) +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(x = "Hazard ratio for waitlist death", y = NULL) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(out_dir, "three_pollutant_hazard_ratios.png"),
  width = 7,
  height = 4,
  dpi = 300
)

message("Wrote primary survival outputs to ", out_dir)
