#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cmprsk)
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

saf_dir <- "/Users/saborpete/Library/CloudStorage/Box-Box/SAF Q2 2025"
pubsaf_dir <- file.path(saf_dir, "pubsaf2506")
supp_dir <- file.path(saf_dir, "SupplementalData2506")
pollution_dir <- "output/zip_pollution"
primary_dir <- "output/primary_survival"
competing_dir <- "output/competing_timevarying"
out_dir <- "output/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analysis_start_year <- 2005L
analysis_end_year <- 2023L
max_curve_years <- 5
curve_times <- seq(0, max_curve_years * 365.25, by = 30.4375)

candidate_files <- tribble(
  ~candidate_group, ~path,
  "kidney_pancreas", file.path(pubsaf_dir, "cand_kipa.sas7bdat"),
  "liver_intestine", file.path(pubsaf_dir, "cand_liin.sas7bdat"),
  "thoracic", file.path(pubsaf_dir, "cand_thor.sas7bdat")
)

clean_zip <- function(x) {
  z <- str_extract(as.character(x), "[0-9]+")
  z <- if_else(nchar(z) > 5L, substr(z, 1, 5), z)
  if_else(is.na(z), NA_character_, str_pad(z, 5, side = "left", pad = "0"))
}

mode_value <- function(x) {
  ux <- na.omit(x)
  if (length(ux) == 0L) return(NA_character_)
  names(which.max(table(ux)))
}

pollutant_labels <- c(
  pm25_ug_m3 = "PM2.5",
  o3_ppb = "Ozone",
  no2 = "NO2"
)
pollutant_order <- c("PM2.5", "Ozone", "NO2")
organ_labels <- c(
  HL = "Heart-lung",
  HR = "Heart",
  IN = "Intestine",
  KI = "Kidney",
  KP = "Kidney-pancreas",
  LI = "Liver",
  LU = "Lung",
  PA = "Pancreas"
)

pollutant_units <- c(
  pm25_ug_m3 = "ug/m3",
  o3_ppb = "ppb",
  no2 = "NO2 units"
)

theme_waitlist <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title.position = "plot",
      legend.position = "bottom"
    )
}

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
  transmute(zip = clean_zip(zip), year, no2 = no2, no2_source = value_source)

pollution_year <- pm25_year %>%
  inner_join(o3_year, by = c("zip", "year")) %>%
  inner_join(no2_year, by = c("zip", "year")) %>%
  filter(
    !is.na(pm25_ug_m3), !is.na(o3_ppb), !is.na(no2),
    pm25_n_months >= 12L, o3_n_months >= 12L
  ) %>%
  mutate(
    pm25_5ug = pm25_ug_m3 / 5,
    o3_10ppb = o3_ppb / 10,
    no2_10unit = no2 / 10
  )

pollution_trends <- pollution_year %>%
  select(year, pm25_ug_m3, o3_ppb, no2) %>%
  tidyr::pivot_longer(-year, names_to = "pollutant", values_to = "value") %>%
  group_by(year, pollutant) %>%
  summarise(
    p25 = quantile(value, 0.25, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    p75 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pollutant = factor(pollutant, levels = names(pollutant_labels), labels = pollutant_labels)
  )

write_csv(pollution_trends, file.path(out_dir, "pollution_trends_by_year.csv"))

p_pollution_trends <- ggplot(pollution_trends, aes(x = year, y = median)) +
  geom_ribbon(aes(ymin = p25, ymax = p75), fill = "#9ecae1", alpha = 0.35) +
  geom_line(color = "#08519c", linewidth = 0.8) +
  facet_wrap(~pollutant, scales = "free_y", nrow = 1) +
  scale_x_continuous(breaks = seq(2005, 2023, by = 3)) +
  labs(
    title = "Annual ZIP-level pollution distributions",
    subtitle = "Median with interquartile range across ZIP-year pollution surfaces",
    x = NULL,
    y = "Annual concentration"
  ) +
  theme_waitlist()

ggsave(file.path(out_dir, "pollution_trends_by_year.png"), p_pollution_trends, width = 10, height = 4.5, dpi = 300)

message("Reading SAF candidate and ZIP files")
candidate_zip <- read_sas(file.path(supp_dir, "canzip2506.sas7bdat")) %>%
  transmute(PERS_ID, PX_ID, candidate_zip = clean_zip(CAN_PERM_ZIP)) %>%
  distinct(PERS_ID, PX_ID, .keep_all = TRUE)

candidate_data <- candidate_files %>%
  rowwise() %>%
  mutate(data = list(
    read_sas(
      path,
      col_select = any_of(c(
        "PERS_ID", "PX_ID", "WL_ORG", "CAN_LISTING_DT", "CAN_ACTIVATE_DT",
        "CAN_SOURCE", "CAN_REM_DT", "CAN_REM_CD", "CAN_DEATH_DT", "CAN_ENDWLFU",
        "CAN_AGE_AT_LISTING", "CAN_GENDER", "CAN_RACE_SRTR", "CAN_ETHNICITY_SRTR",
        "CAN_BMI"
      ))
    ) %>%
      mutate(candidate_group = candidate_group)
  )) %>%
  ungroup() %>%
  select(data) %>%
  tidyr::unnest(data)

message("Constructing completed-analysis cohort")
cohort <- candidate_data %>%
  left_join(candidate_zip, by = c("PERS_ID", "PX_ID")) %>%
  mutate(
    index_date = coalesce(CAN_ACTIVATE_DT, CAN_LISTING_DT),
    index_year = as.integer(format(index_date, "%Y")),
    death_event = CAN_REM_CD == 8,
    transplant_event = CAN_REM_CD %in% c(4, 15, 18, 19, 14, 22, 23),
    event_date = if_else(death_event, coalesce(CAN_DEATH_DT, CAN_REM_DT), CAN_REM_DT),
    censor_date = coalesce(CAN_REM_DT, CAN_ENDWLFU),
    end_date = coalesce(event_date, censor_date),
    followup_days = as.numeric(end_date - index_date),
    event_type = case_when(
      death_event ~ 1L,
      transplant_event ~ 2L,
      TRUE ~ 0L
    ),
    age = CAN_AGE_AT_LISTING,
    sex = factor(if_else(CAN_GENDER %in% c("M", "F"), CAN_GENDER, NA_character_)),
    race_srtr = factor(na_if(CAN_RACE_SRTR, "")),
    ethnicity_srtr = factor(na_if(CAN_ETHNICITY_SRTR, "")),
    wl_org = factor(WL_ORG)
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
  mutate(followup_days = pmax(followup_days, 0.5)) %>%
  left_join(pollution_year, by = c("candidate_zip" = "zip", "index_year" = "year")) %>%
  filter(
    !is.na(pm25_ug_m3), !is.na(o3_ppb), !is.na(no2),
    complete.cases(age, sex, race_srtr, ethnicity_srtr, CAN_BMI, wl_org, index_year)
  ) %>%
  mutate(
    pm25_quartile = ntile(pm25_ug_m3, 4),
    o3_quartile = ntile(o3_ppb, 4),
    no2_quartile = ntile(no2, 4)
  )

cohort_summary <- cohort %>%
  summarise(
    n_registrations = n(),
    n_people = n_distinct(PERS_ID),
    n_deaths = sum(event_type == 1L),
    n_transplants = sum(event_type == 2L),
    n_other_censored = sum(event_type == 0L),
    median_followup_days = median(followup_days)
  )
write_csv(cohort_summary, file.path(out_dir, "figure_cohort_summary.csv"))

make_quartile_table <- function(data, exposure, quartile, label, unit) {
  data %>%
    group_by(quartile = .data[[quartile]]) %>%
    summarise(
      n = n(),
      min_value = min(.data[[exposure]], na.rm = TRUE),
      median_value = median(.data[[exposure]], na.rm = TRUE),
      max_value = max(.data[[exposure]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      pollutant = label,
      unit = unit,
      quartile_label = paste0("Q", quartile, " (", round(min_value, 1), "-", round(max_value, 1), " ", unit, ")")
    )
}

quartile_table <- bind_rows(
  make_quartile_table(cohort, "pm25_ug_m3", "pm25_quartile", "PM2.5", "ug/m3"),
  make_quartile_table(cohort, "o3_ppb", "o3_quartile", "Ozone", "ppb"),
  make_quartile_table(cohort, "no2", "no2_quartile", "NO2", "units")
)
write_csv(quartile_table, file.path(out_dir, "pollution_quartile_ranges.csv"))

extract_km_curve <- function(data, exposure, quartile_var, label) {
  q_lookup <- make_quartile_table(
    data,
    exposure,
    quartile_var,
    label,
    pollutant_units[[exposure]]
  ) %>%
    select(quartile, quartile_label)

  dat <- data %>%
    mutate(exposure_quartile = factor(.data[[quartile_var]], levels = 1:4))

  fit <- survfit(Surv(followup_days, event_type == 1L) ~ exposure_quartile, data = dat)
  s <- summary(fit, times = curve_times, extend = TRUE)

  tibble(
    time_days = s$time,
    years_since_listing = s$time / 365.25,
    stratum = s$strata,
    n_risk = s$n.risk,
    cumulative_death = 1 - s$surv,
    standard_error = s$std.err
  ) %>%
    mutate(
      quartile = as.integer(str_extract(stratum, "[0-9]+")),
      pollutant = factor(label, levels = pollutant_order),
      quartile_group = factor(paste0("Q", quartile), levels = paste0("Q", 1:4))
    ) %>%
    left_join(q_lookup, by = "quartile")
}

km_curves <- bind_rows(
  extract_km_curve(cohort, "pm25_ug_m3", "pm25_quartile", "PM2.5"),
  extract_km_curve(cohort, "o3_ppb", "o3_quartile", "Ozone"),
  extract_km_curve(cohort, "no2", "no2_quartile", "NO2")
)
write_csv(km_curves, file.path(out_dir, "waitlist_death_km_by_pollution_quartile.csv"))

p_km <- ggplot(km_curves, aes(x = years_since_listing, y = cumulative_death, color = quartile_group)) +
  geom_step(linewidth = 0.75) +
  facet_wrap(~pollutant, nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = 0:max_curve_years, limits = c(0, max_curve_years)) +
  labs(
    title = "Cumulative waitlist death by baseline pollution quartile",
    subtitle = "Kaplan-Meier failure functions; transplant and other removal censored",
    x = "Years since waitlist listing/activation",
    y = "Cumulative waitlist death",
    color = "Exposure quartile"
  ) +
  theme_waitlist()

ggsave(file.path(out_dir, "waitlist_death_km_by_pollution_quartile.png"), p_km, width = 12, height = 5, dpi = 300)

extract_cif_curve <- function(data, exposure, quartile_var, label) {
  q_lookup <- make_quartile_table(
    data,
    exposure,
    quartile_var,
    label,
    pollutant_units[[exposure]]
  ) %>%
    select(quartile, quartile_label)

  dat <- data %>%
    mutate(exposure_quartile = .data[[quartile_var]])

  fit <- cuminc(
    ftime = dat$followup_days,
    fstatus = dat$event_type,
    group = dat$exposure_quartile,
    cencode = 0
  )

  bind_rows(lapply(names(fit), function(nm) {
    if (nm == "Tests") return(NULL)
    parts <- str_split(nm, " ", simplify = TRUE)
    quartile <- as.integer(parts[1])
    event_code <- as.integer(parts[2])
    est <- approx(
      x = fit[[nm]]$time,
      y = fit[[nm]]$est,
      xout = curve_times,
      method = "constant",
      rule = 2,
      f = 0
    )$y
    tibble(
      time_days = curve_times,
      years_since_listing = curve_times / 365.25,
      cumulative_incidence = est,
      quartile = quartile,
      event = recode(as.character(event_code), `1` = "Waitlist death", `2` = "Transplant"),
      pollutant = factor(label, levels = pollutant_order),
      quartile_group = factor(paste0("Q", quartile), levels = paste0("Q", 1:4))
    )
  })) %>%
    left_join(q_lookup, by = "quartile")
}

extract_cif_curve_by_organ <- function(data, exposure, quartile_var, label) {
  bind_rows(lapply(sort(unique(as.character(data$WL_ORG))), function(org) {
    organ_data <- data %>% filter(as.character(WL_ORG) == org)
    if (nrow(organ_data) == 0L || length(unique(organ_data[[quartile_var]])) < 2L) return(NULL)

    extract_cif_curve(organ_data, exposure, quartile_var, label) %>%
      mutate(
        WL_ORG = org,
        organ = factor(recode(org, !!!organ_labels), levels = organ_labels)
      )
  }))
}

cif_curves <- bind_rows(
  extract_cif_curve(cohort, "pm25_ug_m3", "pm25_quartile", "PM2.5"),
  extract_cif_curve(cohort, "o3_ppb", "o3_quartile", "Ozone"),
  extract_cif_curve(cohort, "no2", "no2_quartile", "NO2")
)
write_csv(cif_curves, file.path(out_dir, "competing_cumulative_incidence_by_pollution_quartile.csv"))

p_cif_death <- cif_curves %>%
  filter(event == "Waitlist death") %>%
  ggplot(aes(x = years_since_listing, y = cumulative_incidence, color = quartile_group)) +
  geom_step(linewidth = 0.75) +
  facet_wrap(~pollutant, nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = 0:max_curve_years, limits = c(0, max_curve_years)) +
  labs(
    title = "Competing-risk cumulative incidence of waitlist death",
    subtitle = "Transplant treated as a competing event",
    x = "Years since waitlist listing/activation",
    y = "Cumulative incidence",
    color = "Exposure quartile"
  ) +
  theme_waitlist()

ggsave(file.path(out_dir, "waitlist_death_cif_by_pollution_quartile.png"), p_cif_death, width = 12, height = 5, dpi = 300)

organ_cif_curves <- bind_rows(
  extract_cif_curve_by_organ(cohort, "pm25_ug_m3", "pm25_quartile", "PM2.5"),
  extract_cif_curve_by_organ(cohort, "o3_ppb", "o3_quartile", "Ozone"),
  extract_cif_curve_by_organ(cohort, "no2", "no2_quartile", "NO2")
)
write_csv(organ_cif_curves, file.path(out_dir, "organ_competing_cumulative_incidence_by_pollution_quartile.csv"))

for (pollutant_name in pollutant_order) {
  pollutant_slug <- recode(pollutant_name, `PM2.5` = "pm25", Ozone = "o3", NO2 = "no2")

  p_organ_death <- organ_cif_curves %>%
    filter(event == "Waitlist death", pollutant == pollutant_name) %>%
    ggplot(aes(x = years_since_listing, y = cumulative_incidence, color = quartile_group)) +
    geom_step(linewidth = 0.65) +
    facet_wrap(~organ, ncol = 4) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_x_continuous(breaks = 0:max_curve_years, limits = c(0, max_curve_years)) +
    labs(
      title = paste0("Waitlist death cumulative incidence by ", pollutant_name, " quartile"),
      subtitle = "Organ-specific CIFs; transplant treated as a competing event",
      x = "Years since waitlist listing/activation",
      y = "Cumulative incidence",
      color = "Exposure quartile"
    ) +
    theme_waitlist(base_size = 10)

  ggsave(
    file.path(out_dir, paste0("organ_waitlist_death_cif_by_", pollutant_slug, "_quartile.png")),
    p_organ_death,
    width = 12,
    height = 7,
    dpi = 300
  )
}

p_organ_pm25_events <- organ_cif_curves %>%
  filter(pollutant == "PM2.5") %>%
  ggplot(aes(x = years_since_listing, y = cumulative_incidence, color = quartile_group)) +
  geom_step(linewidth = 0.55) +
  facet_grid(event ~ organ, scales = "free_y") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = c(0, 2, 4), limits = c(0, max_curve_years)) +
  labs(
    title = "Organ-specific competing outcomes by PM2.5 quartile",
    subtitle = "Cumulative incidence functions for waitlist death and transplant",
    x = "Years since waitlist listing/activation",
    y = "Cumulative incidence",
    color = "Exposure quartile"
  ) +
  theme_waitlist(base_size = 9) +
  theme(axis.text.x = element_text(size = 7))

ggsave(file.path(out_dir, "organ_pm25_competing_outcomes_cif.png"), p_organ_pm25_events, width = 16, height = 7, dpi = 300)

p_cif_events <- cif_curves %>%
  filter(pollutant == "PM2.5") %>%
  ggplot(aes(x = years_since_listing, y = cumulative_incidence, color = quartile_label)) +
  geom_step(linewidth = 0.75) +
  facet_wrap(~event, nrow = 1, scales = "free_y") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = 0:max_curve_years, limits = c(0, max_curve_years)) +
  labs(
    title = "Competing outcomes by PM2.5 quartile",
    subtitle = "Cumulative incidence functions for completed competing-risk cohort",
    x = "Years since waitlist listing/activation",
    y = "Cumulative incidence",
    color = "Baseline PM2.5"
  ) +
  theme_waitlist()

ggsave(file.path(out_dir, "pm25_competing_outcomes_cif.png"), p_cif_events, width = 9, height = 5, dpi = 300)

message("Making model-summary forest plots")
primary_hr <- read_csv(file.path(primary_dir, "cox_pollutant_results.csv"), show_col_types = FALSE) %>%
  filter(model == "three_pollutant_model") %>%
  transmute(
    model = "Primary Cox: waitlist death",
    term = pollutant,
    estimate = hazard_ratio,
    conf_low,
    conf_high
  )

cause_specific_hr <- read_csv(file.path(competing_dir, "cause_specific_cox_pollutant_results.csv"), show_col_types = FALSE) %>%
  transmute(
    model = recode(
      model,
      cause_specific_waitlist_death = "Cause-specific Cox: waitlist death",
      cause_specific_transplant = "Cause-specific Cox: transplant"
    ),
    term = pollutant,
    estimate = hazard_ratio,
    conf_low,
    conf_high
  )

timevarying_hr <- read_csv(file.path(competing_dir, "timevarying_cox_results.csv"), show_col_types = FALSE) %>%
  transmute(
    model = "Time-varying Cox: waitlist death",
    term,
    estimate = hazard_ratio,
    conf_low,
    conf_high
  )

fine_gray_hr <- read_csv(file.path(competing_dir, "fine_gray_pollutant_results_by_organ.csv"), show_col_types = FALSE) %>%
  filter(!is.na(pollutant)) %>%
  transmute(
    model = paste0(str_replace(model, "fine_gray_waitlist_death_", "Fine-Gray waitlist death: "), " (monitored run)"),
    term = pollutant,
    estimate = subdistribution_hazard_ratio,
    conf_low,
    conf_high
  )

model_hr <- bind_rows(primary_hr, cause_specific_hr, timevarying_hr, fine_gray_hr) %>%
  filter(term %in% c("PM2.5 per 5 ug/m3", "O3 per 10 ppb", "NO2 per 10 units")) %>%
  mutate(
    term = factor(term, levels = rev(c("PM2.5 per 5 ug/m3", "O3 per 10 ppb", "NO2 per 10 units"))),
    model = factor(model, levels = unique(model))
  )
write_csv(model_hr, file.path(out_dir, "completed_model_hazard_ratios_for_plot.csv"))

p_hr <- ggplot(model_hr, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0.18, color = "#525252") +
  geom_point(size = 1.9, color = "#08519c") +
  facet_wrap(~model, scales = "free_x") +
  scale_x_log10() +
  labs(
    title = "Completed model estimates",
    subtitle = "Hazard/subdistribution hazard ratios with 95% confidence intervals",
    x = "Ratio estimate, log scale",
    y = NULL
  ) +
  theme_waitlist(base_size = 10) +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "completed_model_hazard_ratio_forest.png"), p_hr, width = 12, height = 9, dpi = 300)

main_model_hr <- model_hr %>%
  filter(!str_detect(as.character(model), "Fine-Gray")) %>%
  mutate(model = factor(model, levels = c(
    "Primary Cox: waitlist death",
    "Cause-specific Cox: waitlist death",
    "Cause-specific Cox: transplant",
    "Time-varying Cox: waitlist death"
  )))

p_main_hr <- ggplot(main_model_hr, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0.18, color = "#525252") +
  geom_point(size = 2.2, color = "#08519c") +
  facet_wrap(~model, scales = "free_x", nrow = 1) +
  scale_x_log10() +
  labs(
    title = "Main completed model estimates",
    subtitle = "Adjusted hazard ratios with 95% confidence intervals",
    x = "Hazard ratio, log scale",
    y = NULL
  ) +
  theme_waitlist(base_size = 11) +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "main_completed_model_hazard_ratio_forest.png"), p_main_hr, width = 12, height = 4.5, dpi = 300)

fine_gray_model_hr <- model_hr %>%
  filter(str_detect(as.character(model), "Fine-Gray"))

p_fg_hr <- ggplot(fine_gray_model_hr, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0.18, color = "#525252") +
  geom_point(size = 1.9, color = "#08519c") +
  facet_wrap(~model, scales = "free_x", ncol = 4) +
  scale_x_log10() +
  labs(
    title = "Organ-specific Fine-Gray estimates from monitored run",
    subtitle = "Subdistribution hazard ratios with 95% confidence intervals",
    x = "Subdistribution hazard ratio, log scale",
    y = NULL
  ) +
  theme_waitlist(base_size = 10) +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "fine_gray_monitored_hazard_ratio_forest.png"), p_fg_hr, width = 12, height = 7, dpi = 300)

message("Wrote figures to ", out_dir)
