#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

source(file.path("code", "rolling_prior_pollution.R"))

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(sf)
  library(stringr)
  library(survival)
  library(tibble)
})

sf_use_s2(FALSE)

analysis_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_analysis_dataset.csv.gz"
)
community_path <- file.path("data", "processed", "community", "zcta_acs_community_covariates_2005_2023.csv.gz")
release_dir <- file.path("data", "release")
pollution_dir <- file.path(release_dir, "air_pollution_zcta_parquet")
zcta_path <- file.path("data", "cache", "cb_2020_us_zcta520_500k", "cb_2020_us_zcta520_500k.shp")
prior_cache_dir <- file.path("output", "prior_year_pollution_cox_svi", "cache")
fig_dir <- file.path("output", "figures", "pollution_absolute_1yr_survival_maps")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(analysis_path)) stop("Missing analysis dataset: ", analysis_path, call. = FALSE)
if (!file.exists(zcta_path)) stop("Missing cached ZCTA shapefile: ", zcta_path, call. = FALSE)

organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")
target_organs <- names(organ_labels)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

clean_zip <- function(x) {
  z <- str_extract(as.character(x), "[0-9]+")
  z <- if_else(nchar(z) > 5L, substr(z, 1, 5), z)
  if_else(is.na(z), NA_character_, str_pad(z, 5, side = "left", pad = "0"))
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  names(sort(table(as.character(x)), decreasing = TRUE))[1]
}

make_complete_acs_svi_proxy <- function(path) {
  community <- read_csv(path, show_col_types = FALSE) %>%
    mutate(
      zip = clean_zip(zip),
      analysis_year = as.integer(analysis_year)
    )

  community_2023 <- community %>%
    filter(analysis_year == 2023L) %>%
    select(-analysis_year)

  community <- bind_rows(
    community,
    community_2023 %>% mutate(analysis_year = 2024L),
    community_2023 %>% mutate(analysis_year = 2025L)
  )

  vulnerability_vars <- c(
    "pct_poverty", "pct_unemployed", "pct_no_vehicle", "pct_nonwhite",
    "median_household_income", "pct_bachelor_plus"
  )

  community %>%
    group_by(analysis_year) %>%
    mutate(across(
      all_of(vulnerability_vars),
      ~ if_else(is.na(.x), median(.x, na.rm = TRUE), .x)
    )) %>%
    mutate(
      svi_poverty_rank = percent_rank(pct_poverty),
      svi_unemployed_rank = percent_rank(pct_unemployed),
      svi_no_vehicle_rank = percent_rank(pct_no_vehicle),
      svi_nonwhite_rank = percent_rank(pct_nonwhite),
      svi_low_income_rank = percent_rank(-median_household_income),
      svi_low_education_rank = percent_rank(-pct_bachelor_plus),
      zcta_svi_proxy = rowMeans(
        cbind(
          svi_poverty_rank,
          svi_unemployed_rank,
          svi_no_vehicle_rank,
          svi_nonwhite_rank,
          svi_low_income_rank,
          svi_low_education_rank
        ),
        na.rm = TRUE
      )
    ) %>%
    ungroup() %>%
    select(zip, analysis_year, zcta_svi_proxy)
}

read_map_annual_pollution <- function() {
  pm25 <- read_daily_pollution_aggregate(
    file.path(release_dir, "lghap_pm25_zcta_daily_parquet"),
    "pm25_ug_m3", "pm25_prior_ug_m3", "annual",
    file.path(prior_cache_dir, "pm25_daily_annual_zcta.csv.gz")
  )

  no2 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_no2_annual_2005_2025.parquet")) %>%
    transmute(zip = clean_zip(zip), year = as.integer(year), no2_prior_ppb = no2)

  list(pm25 = pm25, no2 = no2)
}

read_pollutant_mean <- function(name, spec) {
  if (name == "pm25") {
    return(read_map_annual_pollution()$pm25 %>%
      filter(year >= spec$start_year, year <= spec$end_year) %>%
      group_by(zip) %>%
      summarise(value = mean(pm25_prior_ug_m3), years = n(), .groups = "drop") %>%
      filter(years == (spec$end_year - spec$start_year + 1L)) %>%
      mutate(exposure = name))
  }
  dat <- read_parquet(file.path(pollution_dir, spec$parquet)) %>%
    transmute(
      zip = clean_zip(zip),
      year = as.integer(year),
      value = .data[[spec$value_col]]
    ) %>%
    filter(year >= spec$start_year, year <= spec$end_year)

  if (identical(spec$table, "monthly")) {
    dat %>%
      group_by(zip, year) %>%
      summarise(n_periods = sum(!is.na(value)), annual_mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
      filter(n_periods == 12L, is.finite(annual_mean)) %>%
      group_by(zip) %>%
      summarise(value = mean(annual_mean, na.rm = TRUE), years = n(), .groups = "drop") %>%
      filter(years == (spec$end_year - spec$start_year + 1L)) %>%
      mutate(exposure = name)
  } else {
    dat %>%
      filter(is.finite(value)) %>%
      group_by(zip) %>%
      summarise(value = mean(value, na.rm = TRUE), years = n(), .groups = "drop") %>%
      filter(years == (spec$end_year - spec$start_year + 1L)) %>%
      mutate(exposure = name)
  }
}

baseline_hazard_at <- function(fit, model_dat, time_days = 365) {
  bh <- basehaz(fit, centered = FALSE)
  center_weights <- model_dat %>%
    count(listing_center, name = "n") %>%
    mutate(
      listing_center = as.character(listing_center),
      center_weight = n / sum(n)
    )

  if (!"strata" %in% names(bh)) {
    return(max(bh$hazard[bh$time <= time_days], na.rm = TRUE))
  }

  bh_at <- bh %>%
    filter(time <= time_days) %>%
    group_by(strata) %>%
    slice_max(time, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      listing_center = sub("^.*=", "", as.character(strata)),
      hazard_365 = hazard
    )

  center_weights %>%
    left_join(bh_at, by = "listing_center") %>%
    mutate(hazard_365 = coalesce(hazard_365, 0)) %>%
    summarise(hazard_365 = sum(center_weight * hazard_365), .groups = "drop") %>%
    pull(hazard_365)
}

fit_survival_map_model <- function(dat, zcta_pollution, org, pollutant, term, unit_label) {
  vars_needed <- c("followup_days", "adverse_event", term, "age", "sex", "race", "zcta_svi_proxy", "listing_center")
  model_dat <- dat %>%
    filter(
      WL_ORG == org,
      complete.cases(across(all_of(vars_needed))),
      followup_days > 0
    ) %>%
    mutate(
      .followup_days = pmax(followup_days, 0.5),
      sex = droplevels(factor(sex)),
      race = droplevels(factor(race)),
      listing_center = droplevels(factor(listing_center))
    )

  log_msg(
    "Fitting absolute-scale Cox model ",
    recode(org, !!!organ_labels),
    " | ",
    pollutant,
    " n=",
    nrow(model_dat),
    " adverse=",
    sum(model_dat$adverse_event)
  )

  form <- as.formula(paste(
    "Surv(.followup_days, adverse_event) ~",
    paste(c(term, "age", "sex", "race", "zcta_svi_proxy", "strata(listing_center)"), collapse = " + ")
  ))
  fit <- coxph(form, data = model_dat, ties = "efron", x = TRUE, y = TRUE)
  h0_365 <- baseline_hazard_at(fit, model_dat, time_days = 365)

  ref <- tibble(
    age = median(model_dat$age, na.rm = TRUE),
    sex = factor(mode_value(model_dat$sex), levels = levels(model_dat$sex)),
    race = factor(mode_value(model_dat$race), levels = levels(model_dat$race)),
    zcta_svi_proxy = median(model_dat$zcta_svi_proxy, na.rm = TRUE),
    listing_center = factor(levels(model_dat$listing_center)[1], levels = levels(model_dat$listing_center))
  )

  newdat <- zcta_pollution %>%
    filter(exposure == pollutant, is.finite(value)) %>%
    transmute(
      zip,
      exposure = pollutant,
      value,
      !!term := value / if_else(pollutant == "pm25", 5, 10),
      age = ref$age,
      sex = factor(as.character(ref$sex), levels = levels(model_dat$sex)),
      race = factor(as.character(ref$race), levels = levels(model_dat$race)),
      zcta_svi_proxy = ref$zcta_svi_proxy,
      listing_center = factor(as.character(ref$listing_center), levels = levels(model_dat$listing_center))
    )

  lp <- predict(fit, newdata = newdat, type = "lp", reference = "zero")
  survival_1yr <- exp(-h0_365 * exp(lp))

  tibble(
    zip = newdat$zip,
    organ = recode(org, !!!organ_labels),
    exposure = pollutant,
    exposure_value = newdat$value,
    survival_1yr = survival_1yr,
    adverse_risk_1yr = 1 - survival_1yr,
    baseline_hazard_365 = h0_365,
    reference_age = ref$age,
    reference_sex = as.character(ref$sex),
    reference_race = as.character(ref$race),
    reference_svi = ref$zcta_svi_proxy,
    model_n = nrow(model_dat),
    adverse_events = sum(model_dat$adverse_event),
    exposure_unit = unit_label
  )
}

pollutant_specs <- list(
  pm25 = list(
    label = "PM2.5",
    unit_label = "per 5 \u00b5g/m\u00b3",
    term = "pm25_prior_5ug",
    value_col = "pm25_ug_m3",
    table = "daily",
    start_year = 2005L,
    end_year = 2024L
  ),
  no2 = list(
    label = "NO2",
    unit_label = "per 10 ppb",
    term = "no2_prior_10ppb",
    parquet = "air_pollution_zcta_no2_annual_2005_2025.parquet",
    value_col = "no2",
    table = "annual",
    start_year = 2006L,
    end_year = 2025L
  )
)

pollutant_panel_labels <- c(pm25 = "PM[2.5]", no2 = "NO[2]")
organ_panel_labels <- c(Kidney = "Kidney", Liver = "Liver", Heart = "Heart", Lung = "Lung")
risk_palette <- c(
  "#1B1B3A", "#243B6B", "#1F5F8B", "#1F7A8C", "#2A9D8F",
  "#7BC8A4", "#B7E4A8", "#E9F5A1"
)

log_msg("Reading primary deduplicated cohort")
analysis_dat <- read_csv(analysis_path, show_col_types = FALSE) %>%
  mutate(
    candidate_zip = clean_zip(candidate_zip),
    index_date = as.Date(index_date),
    observed_end_date = as.Date(observed_end_date),
    listing_year_int = as.integer(as.character(listing_year)),
    followup_days = as.numeric(observed_end_date - index_date),
    sex = factor(sex),
    race = factor(race),
    listing_center = factor(listing_center)
  ) %>%
  filter(WL_ORG %in% target_organs)

log_msg("Attaching ACS-derived ZCTA SVI proxy")
svi <- make_complete_acs_svi_proxy(community_path)
analysis_dat <- analysis_dat %>%
  left_join(svi, by = c("candidate_zip" = "zip", "listing_year_int" = "analysis_year"), suffix = c("", "_community")) %>%
  mutate(zcta_svi_proxy = coalesce(zcta_svi_proxy, zcta_svi_proxy_community)) %>%
  select(-any_of("zcta_svi_proxy_community"))

log_msg("Attaching prior-year pollution")
prior_pollution <- read_rolling_prior_pollution(analysis_dat)
analysis_dat <- analysis_dat %>%
  left_join(prior_pollution$pm25, by = c("candidate_zip" = "zip", "index_date" = "index_date")) %>%
  left_join(prior_pollution$no2, by = c("candidate_zip" = "zip", "index_date" = "index_date")) %>%
  mutate(
    pm25_prior_5ug = pm25_prior_ug_m3 / 5,
    no2_prior_10ppb = no2_prior_ppb / 10
  )

log_msg("Reading study-period ZCTA pollution means")
pollution_means <- bind_rows(lapply(names(pollutant_specs), function(name) {
  read_pollutant_mean(name, pollutant_specs[[name]])
}))

map_values <- bind_rows(lapply(target_organs, function(org) {
  bind_rows(lapply(names(pollutant_specs), function(pollutant) {
    spec <- pollutant_specs[[pollutant]]
    fit_survival_map_model(
      analysis_dat,
      pollution_means,
      org = org,
      pollutant = pollutant,
      term = spec$term,
      unit_label = spec$unit_label
    )
  }))
}))

write_csv(map_values, file.path(fig_dir, "pollution_absolute_1yr_survival_map_values.csv"))

summary_tbl <- map_values %>%
  group_by(organ, exposure) %>%
  summarise(
    zctas = n(),
    model_n = first(model_n),
    adverse_events = first(adverse_events),
    baseline_hazard_365 = first(baseline_hazard_365),
    reference_age = first(reference_age),
    reference_sex = first(reference_sex),
    reference_race = first(reference_race),
    reference_svi = first(reference_svi),
    min_survival_1yr = min(survival_1yr, na.rm = TRUE),
    p01_survival_1yr = as.numeric(quantile(survival_1yr, 0.01, na.rm = TRUE)),
    median_survival_1yr = median(survival_1yr, na.rm = TRUE),
    p99_survival_1yr = as.numeric(quantile(survival_1yr, 0.99, na.rm = TRUE)),
    max_survival_1yr = max(survival_1yr, na.rm = TRUE),
    min_adverse_risk_1yr = min(adverse_risk_1yr, na.rm = TRUE),
    p99_adverse_risk_1yr = as.numeric(quantile(adverse_risk_1yr, 0.99, na.rm = TRUE)),
    max_adverse_risk_1yr = max(adverse_risk_1yr, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(summary_tbl, file.path(fig_dir, "pollution_absolute_1yr_survival_map_summary.csv"))

log_msg("Reading CONUS ZCTA boundaries")
zcta <- st_read(zcta_path, quiet = TRUE)[, c("ZCTA5CE20", "geometry")]
names(zcta)[1] <- "zip"
zcta$zip <- clean_zip(zcta$zip)
conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326))
zcta <- suppressWarnings(st_crop(zcta, st_transform(conus_bbox, st_crs(zcta))))
zcta <- st_transform(zcta, 5070)

plot_dat <- merge(
  zcta,
  map_values %>%
    transmute(
      zip,
      exposure,
      pollutant_panel = factor(exposure, levels = names(pollutant_specs)),
      organ_panel = factor(organ_panel_labels[organ], levels = organ_panel_labels),
      adverse_risk_1yr_pct = adverse_risk_1yr * 100
    ),
  by = "zip",
  all.x = FALSE
)

scale_limits <- plot_dat %>%
  st_drop_geometry() %>%
  summarise(
    lo = floor(as.numeric(quantile(adverse_risk_1yr_pct, 0.01, na.rm = TRUE))),
    hi = ceiling(as.numeric(quantile(adverse_risk_1yr_pct, 0.99, na.rm = TRUE))),
    .groups = "drop"
  )
legend_breaks <- pretty(c(scale_limits$lo, scale_limits$hi), n = 8)

plot_dat <- plot_dat %>%
  group_by(organ_panel) %>%
  mutate(
    organ_scale_lo = floor(as.numeric(quantile(adverse_risk_1yr_pct, 0.01, na.rm = TRUE)) * 10) / 10,
    organ_scale_hi = ceiling(as.numeric(quantile(adverse_risk_1yr_pct, 0.99, na.rm = TRUE)) * 10) / 10,
    organ_risk_winsor = pmin(pmax(adverse_risk_1yr_pct, organ_scale_lo), organ_scale_hi),
    organ_risk_scaled = scales::rescale(organ_risk_winsor, to = c(0, 1), from = c(first(organ_scale_lo), first(organ_scale_hi)))
  ) %>%
  ungroup()

write_csv(
  tibble(
    scale = "shared_1yr_adverse_event_risk_percent_winsorized_at_global_p01_p99",
    lo = scale_limits$lo,
    hi = scale_limits$hi,
    breaks = paste(sprintf("%.0f", legend_breaks), collapse = ", ")
  ),
  file.path(fig_dir, "pollution_absolute_1yr_survival_shared_scale.csv")
)

theme_map <- function() {
  theme_void(base_size = 13) +
    theme(
      plot.title = element_blank(),
      strip.text = element_text(face = "bold", size = 16, color = "grey8", margin = margin(8, 8, 8, 8)),
      strip.text.x = element_text(face = "bold"),
      strip.text.y.right = element_text(angle = 270),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.55),
      legend.position = "bottom",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.key.width = unit(160, "pt"),
      legend.key.height = unit(16, "pt"),
      panel.spacing = unit(5, "pt"),
      plot.margin = margin(5, 8, 5, 8)
    )
}

panel <- ggplot(plot_dat) +
  geom_sf(aes(fill = organ_risk_scaled), color = NA) +
  facet_grid(
    organ_panel ~ pollutant_panel,
    labeller = labeller(pollutant_panel = as_labeller(pollutant_panel_labels, label_parsed))
  ) +
  scale_fill_gradientn(
    colours = risk_palette,
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("Lower", "25%", "50%", "75%", "Higher"),
    oob = scales::squish,
    name = "Within-organ relative level of model-implied 1-year risk",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(16, "pt"),
      barwidth = unit(580, "pt"),
      ticks = TRUE,
      ticks.linewidth = 0.35
    )
  ) +
  coord_sf(datum = NA) +
  theme_map()

png_path <- file.path(fig_dir, "figure_pollution_absolute_1yr_survival_maps.png")
pdf_path <- file.path(fig_dir, "figure_pollution_absolute_1yr_survival_maps.pdf")

ggsave(png_path, panel, width = 10.5, height = 13.5, dpi = 360, device = grDevices::png, bg = "white")
ggsave(pdf_path, panel, width = 10.5, height = 13.5, device = grDevices::pdf, bg = "white")

make_organ_panel <- function(organ_name) {
  organ_dat <- plot_dat %>%
    filter(as.character(organ_panel) == organ_name)

  organ_limits <- organ_dat %>%
    st_drop_geometry() %>%
    summarise(
      lo = floor(as.numeric(quantile(adverse_risk_1yr_pct, 0.01, na.rm = TRUE)) * 10) / 10,
      hi = ceiling(as.numeric(quantile(adverse_risk_1yr_pct, 0.99, na.rm = TRUE)) * 10) / 10,
      .groups = "drop"
    )
  organ_breaks <- pretty(c(organ_limits$lo, organ_limits$hi), n = 7)

  organ_dat <- organ_dat %>%
    mutate(organ_risk_winsor = pmin(pmax(adverse_risk_1yr_pct, organ_limits$lo), organ_limits$hi))

  organ_plot <- ggplot(organ_dat) +
    geom_sf(aes(fill = organ_risk_winsor), color = NA) +
    facet_wrap(
      vars(pollutant_panel),
      nrow = 1,
      labeller = labeller(pollutant_panel = as_labeller(pollutant_panel_labels, label_parsed))
    ) +
    scale_fill_gradientn(
      colours = risk_palette,
      limits = c(organ_limits$lo, organ_limits$hi),
      breaks = organ_breaks,
      labels = function(x) sprintf("%.1f%%", x),
      oob = scales::squish,
      name = paste0(organ_name, " model-implied 1-year risk"),
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        barheight = unit(16, "pt"),
        barwidth = unit(520, "pt"),
        ticks = TRUE,
        ticks.linewidth = 0.35
      )
    ) +
    coord_sf(datum = NA) +
    theme_void(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold", size = 18, color = "grey8", margin = margin(8, 8, 8, 8)),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.55),
      legend.position = "bottom",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      panel.spacing = unit(6, "pt"),
      plot.margin = margin(8, 8, 8, 8)
    )

  slug <- tolower(organ_name)
  organ_png <- file.path(fig_dir, paste0("figure_pollution_absolute_1yr_risk_maps_", slug, ".png"))
  organ_pdf <- file.path(fig_dir, paste0("figure_pollution_absolute_1yr_risk_maps_", slug, ".pdf"))
  ggsave(organ_png, organ_plot, width = 10.5, height = 6.2, dpi = 360, device = grDevices::png, bg = "white")
  ggsave(organ_pdf, organ_plot, width = 10.5, height = 6.2, device = grDevices::pdf, bg = "white")

  tibble(
    organ = organ_name,
    png_path = organ_png,
    pdf_path = organ_pdf,
    scale_lo = organ_limits$lo,
    scale_hi = organ_limits$hi,
    scale_breaks = paste(sprintf("%.1f", organ_breaks), collapse = ", ")
  )
}

organ_figure_manifest <- bind_rows(lapply(names(organ_panel_labels), make_organ_panel))
write_csv(organ_figure_manifest, file.path(fig_dir, "pollution_absolute_1yr_risk_organ_figure_manifest.csv"))

write_csv(
  tibble(
    figure_path = c(png_path, pdf_path, organ_figure_manifest$png_path, organ_figure_manifest$pdf_path),
    description = paste(
      "ZCTA maps of model-implied 1-year risk of death/deterioration.",
      "The all-organ panel displays within-organ relative risk levels to emphasize geographic pattern despite different baseline failure rates across organs; individual organ panels retain absolute percent-risk legends.",
      "Predictions use single-pollutant prior-year cause-specific Cox models adjusted for age, sex, race, ZCTA SVI proxy, and center-stratified baseline hazards.",
      "Covariates are fixed at organ-specific reference values and the 365-day baseline cumulative hazard is averaged across listing-center strata by model-cohort center distribution."
    )
  ),
  file.path(fig_dir, "pollution_absolute_1yr_survival_map_manifest.csv")
)

log_msg("Wrote absolute-scale 1-year survival maps to ", normalizePath(fig_dir, winslash = "/"))
