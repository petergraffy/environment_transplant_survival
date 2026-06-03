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
  library(readr)
  library(scales)
  library(stringr)
  library(tidyr)
})

out_dir <- file.path("output", "diagnostics", "pollution_attenuation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pollution_dir <- file.path("output", "zip_pollution")
community_dir <- file.path("data", "processed", "community")
figure_dir <- file.path("output", "figures", "adjusted_competing_risks")

clean_zip <- function(x) {
  z <- str_extract(as.character(x), "[0-9]+")
  z <- if_else(nchar(z) > 5L, substr(z, 1, 5), z)
  if_else(is.na(z), NA_character_, str_pad(z, 5, side = "left", pad = "0"))
}

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
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
    filter(pm25_n_months >= 12L, o3_n_months >= 12L, !is.na(pm25_ug_m3), !is.na(o3_ppb), !is.na(no2))
}

tidy_r2 <- function(model_name, fit) {
  tibble(model = model_name, r_squared = summary(fit)$r.squared, adj_r_squared = summary(fit)$adj.r.squared)
}

log_msg("Reading ZIP-year exposure and community data")
pollution <- read_pollution_annual() %>%
  filter(year >= 2005L, year <= 2023L) %>%
  mutate(zip_prefix = factor(substr(zip, 1, 1)), year_factor = factor(year))

community <- read_csv(
  file.path(community_dir, "zcta_acs_community_covariates_2005_2023.csv.gz"),
  col_types = cols(zip = col_character(), .default = col_guess()),
  show_col_types = FALSE
) %>%
  transmute(
    zip = clean_zip(zip),
    year = analysis_year,
    median_household_income,
    pct_poverty,
    pct_bachelor_plus,
    pct_unemployed,
    pct_no_vehicle,
    pct_black,
    pct_hispanic
  )

zip_year <- pollution %>%
  left_join(community, by = c("zip", "year")) %>%
  filter(complete.cases(median_household_income, pct_poverty, pct_bachelor_plus, pct_unemployed, pct_no_vehicle, pct_black, pct_hispanic))

pollutants <- tribble(
  ~pollutant, ~value_col, ~unit,
  "PM2.5", "pm25_ug_m3", "ug/m3",
  "Ozone", "o3_ppb", "ppb",
  "NO2", "no2", "units"
)

community_terms <- paste(
  "median_household_income + pct_poverty + pct_bachelor_plus +",
  "pct_unemployed + pct_no_vehicle + pct_black + pct_hispanic"
)

log_msg("Estimating how much ZIP-year exposure variation is explained by time, geography, and community covariates")
exposure_r2 <- lapply(seq_len(nrow(pollutants)), function(i) {
  value_col <- pollutants$value_col[[i]]
  pollutant <- pollutants$pollutant[[i]]
  dat <- zip_year %>% filter(!is.na(.data[[value_col]]))
  bind_rows(
    tidy_r2("calendar_year", lm(reformulate("year_factor", value_col), data = dat)),
    tidy_r2("zip_prefix", lm(reformulate("zip_prefix", value_col), data = dat)),
    tidy_r2("year_plus_zip_prefix", lm(reformulate(c("year_factor", "zip_prefix"), value_col), data = dat)),
    tidy_r2("community", lm(as.formula(paste(value_col, "~", community_terms)), data = dat)),
    tidy_r2("year_plus_community", lm(as.formula(paste(value_col, "~ year_factor +", community_terms)), data = dat)),
    tidy_r2("year_plus_zip_prefix_plus_community", lm(as.formula(paste(value_col, "~ year_factor + zip_prefix +", community_terms)), data = dat))
  ) %>%
    mutate(pollutant = pollutant, n_zip_years = nrow(dat), .before = 1)
}) %>%
  bind_rows()

write_csv(exposure_r2, file.path(out_dir, "exposure_explained_variance_zip_year.csv"))

pollution_correlations <- zip_year %>%
  summarise(
    pm25_o3 = cor(pm25_ug_m3, o3_ppb, use = "complete.obs"),
    pm25_no2 = cor(pm25_ug_m3, no2, use = "complete.obs"),
    o3_no2 = cor(o3_ppb, no2, use = "complete.obs")
  ) %>%
  pivot_longer(everything(), names_to = "pollutant_pair", values_to = "pearson_r")

write_csv(pollution_correlations, file.path(out_dir, "pollutant_pairwise_correlations_zip_year.csv"))

log_msg("Summarising 5-year adjusted CIF contrasts from current figures")
curve_files <- c(
  PM2.5 = file.path(figure_dir, "adjusted_competing_risk_curves_pm25.csv"),
  Ozone = file.path(figure_dir, "adjusted_competing_risk_curves_o3.csv"),
  NO2 = file.path(figure_dir, "adjusted_competing_risk_curves_no2.csv")
)

five_year_contrasts <- lapply(names(curve_files), function(pollutant) {
  read_csv(curve_files[[pollutant]], show_col_types = FALSE) %>%
    group_by(organ, event, quartile) %>%
    slice_max(years_since_listing, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(organ, event, quartile, cumulative_incidence) %>%
    pivot_wider(names_from = quartile, values_from = cumulative_incidence) %>%
    mutate(
      pollutant = pollutant,
      q4_minus_q1 = Q4 - Q1,
      relative_q4_vs_q1 = Q4 / Q1,
      .before = 1
    )
}) %>%
  bind_rows()

write_csv(five_year_contrasts, file.path(out_dir, "adjusted_cif_5year_q4_vs_q1_contrasts.csv"))

log_msg("Comparing pollutant estimates before and after benchmark/community adjustment")
primary_path <- file.path("output", "primary_survival", "cox_pollutant_results_by_organ.csv")
score_path <- file.path("output", "score_benchmark_pollution", "score_benchmark_pollution_coefficients.csv")

primary <- read_csv(primary_path, show_col_types = FALSE) %>%
  filter(WL_ORG %in% c("HR", "KI", "LI", "LU")) %>%
  transmute(
    organ = recode(WL_ORG, HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung"),
    pollutant = recode(term, pm25_5ug = "PM2.5", o3_10ppb = "Ozone", no2_10unit = "NO2"),
    model_stage = "original_primary",
    estimate,
    conf.low,
    conf.high,
    p.value
  )

score_prior <- read_csv(score_path, show_col_types = FALSE) %>%
  filter(
    organ %in% c("HR", "KI", "LI", "LU"),
    model %in% c(
      "score_prior_pollution",
      "score_year_prior_pollution",
      "score_year_community_prior_pollution",
      "score_year_center_prior_pollution",
      "score_center_year_community_prior_pollution"
    ),
    str_detect(term, "prior_1y")
  ) %>%
  transmute(
    organ = recode(organ, HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung"),
    pollutant = case_when(
      str_detect(term, "pm25") ~ "PM2.5",
      str_detect(term, "o3") ~ "Ozone",
      str_detect(term, "no2") ~ "NO2",
      TRUE ~ term
    ),
    model_stage = model,
    estimate,
    conf.low,
    conf.high,
    p.value
  )

score_current <- read_csv(score_path, show_col_types = FALSE) %>%
  filter(
    organ %in% c("HR", "KI", "LI", "LU"),
    model == "score_center_year_community_prior_current_pollution",
    str_detect(term, "current")
  ) %>%
  transmute(
    organ = recode(organ, HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung"),
    pollutant = case_when(
      str_detect(term, "pm25") ~ "PM2.5",
      str_detect(term, "o3") ~ "Ozone",
      str_detect(term, "no2") ~ "NO2",
      TRUE ~ term
    ),
    model_stage = "score_year_community_current",
    estimate,
    conf.low,
    conf.high,
    p.value
  )

attenuation <- bind_rows(primary, score_prior, score_current) %>%
  arrange(organ, pollutant, factor(model_stage, levels = c(
    "original_primary",
    "score_prior_pollution",
    "score_year_prior_pollution",
    "score_year_community_prior_pollution",
    "score_year_center_prior_pollution",
    "score_center_year_community_prior_pollution",
    "score_year_community_current"
  )))

write_csv(attenuation, file.path(out_dir, "pollutant_estimate_attenuation_summary.csv"))

attenuation_plot <- attenuation %>%
  mutate(model_stage = factor(model_stage, levels = c(
    "original_primary",
    "score_prior_pollution",
    "score_year_prior_pollution",
    "score_year_community_prior_pollution",
    "score_year_center_prior_pollution",
    "score_center_year_community_prior_pollution",
    "score_year_community_current"
  )))

p1 <- ggplot(attenuation_plot, aes(x = model_stage, y = estimate, ymin = conf.low, ymax = conf.high, color = pollutant)) +
  geom_hline(yintercept = 1, linewidth = 0.35, linetype = "dashed", color = "grey45") +
  geom_pointrange(position = position_dodge(width = 0.45), linewidth = 0.45) +
  facet_wrap(~organ, ncol = 2) +
  scale_y_log10() +
  labs(x = NULL, y = "Hazard ratio / submodel estimate", color = "Pollutant") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "pollutant_estimate_attenuation_summary.png"), p1, width = 10, height = 7, dpi = 300)

p2 <- exposure_r2 %>%
  filter(model %in% c("calendar_year", "zip_prefix", "year_plus_community", "year_plus_zip_prefix_plus_community")) %>%
  ggplot(aes(x = model, y = adj_r_squared, fill = pollutant)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Adjusted R-squared", fill = "Pollutant") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "exposure_explained_variance_zip_year.png"), p2, width = 9, height = 5.5, dpi = 300)

log_msg("Wrote diagnostics to ", out_dir)
