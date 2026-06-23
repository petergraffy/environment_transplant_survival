#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(cmprsk)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

in_path <- file.path("output", "formal_waitlist_environment", "formal_waitlist_environment_analysis_dataset.csv.gz")
out_dir <- file.path("output", "formal_waitlist_environment_finegray")
fig_dir <- file.path("output", "figures", "formal_waitlist_environment_finegray")
model_dir <- file.path(out_dir, "model_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

target_organs <- c("HR", "KI", "LI", "LU")
organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")
finegray_max_n <- as.integer(Sys.getenv("FINEGRAY_ALL_MAX_N", "10000"))

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA)
  names(which.max(table(as.character(x))))
}

collapse_top_levels <- function(x, n_top = 20L) {
  y <- as.character(x)
  y[is.na(y) | y == ""] <- "Missing"
  top <- names(sort(table(y), decreasing = TRUE))[seq_len(min(n_top, length(unique(y))))]
  factor(if_else(y %in% top, y, "Other"))
}

theme_cif <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

model_specs <- list(
  heat_humidity = c("tmax_5c", "rmax_10pct"),
  pm25 = "pm25_5ug",
  o3 = "o3_10ppb",
  no2 = "no2_10unit"
)

curve_terms <- tribble(
  ~model, ~term, ~exposure_label,
  "heat_humidity", "tmax_5c", "Tmax",
  "heat_humidity", "rmax_10pct", "Maximum relative humidity",
  "pm25", "pm25_5ug", "PM2.5",
  "o3", "o3_10ppb", "Ozone",
  "no2", "no2_10unit", "NO2"
)

term_labels <- c(
  tmax_5c = "Tmax per 5 C",
  rmax_10pct = "Maximum relative humidity per 10 pct",
  pm25_5ug = "PM2.5 per 5 ug/m3",
  o3_10ppb = "Ozone per 10 ppb",
  no2_10unit = "NO2 per 10 units"
)

log_msg("Reading formal static analysis dataset")
base <- read_csv(in_path, show_col_types = FALSE) %>%
  filter(WL_ORG %in% target_organs, followup_days > 0) %>%
  mutate(
    fstatus = case_when(
      event_type == "adverse" ~ 1L,
      event_type %in% c("transplant_or_improvement", "other_exit") ~ 2L,
      TRUE ~ 0L
    ),
    sex = factor(sex),
    race = factor(race),
    listing_center = factor(listing_center),
    tmax_5c = tmax_c_waitlist_mean_exact_daily / 5,
    rmax_10pct = rmax_pct_waitlist_mean_exact_daily / 10,
    pm25_5ug = pm25_waitlist_ug_m3 / 5,
    o3_10ppb = o3_waitlist_ppb / 10,
    no2_10unit = no2_waitlist / 10
  )

make_reference_rows <- function(model_dat, exposure_terms, focal_term) {
  ref <- tibble(
    age = median(model_dat$age, na.rm = TRUE),
    sex = factor(mode_value(model_dat$sex), levels = levels(model_dat$sex)),
    race = factor(mode_value(model_dat$race), levels = levels(model_dat$race)),
    index_year_centered = median(model_dat$index_year_centered, na.rm = TRUE),
    organ_score = median(model_dat$organ_score, na.rm = TRUE),
    listing_center = factor(mode_value(model_dat$listing_center), levels = levels(model_dat$listing_center))
  )
  for (term in exposure_terms) {
    ref[[term]] <- median(model_dat[[term]], na.rm = TRUE)
  }
  qs <- quantile(model_dat[[focal_term]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  bind_rows(lapply(seq_along(qs), function(i) {
    out <- ref
    out[[focal_term]] <- qs[[i]]
    out$curve_group <- paste0("Q", c(1, 2, 4)[[i]], " reference")
    out$curve_exposure_value <- qs[[i]]
    out
  }))
}

fit_finegray_model <- function(dat, org, model_name, exposure_terms) {
  vars_needed <- c("followup_days", "fstatus", "age", "sex", "race", "index_year_centered", "organ_score", "listing_center", exposure_terms)
  source_dat <- dat %>%
    filter(complete.cases(across(all_of(vars_needed)))) %>%
    droplevels()
  if (nrow(source_dat) == 0L || sum(source_dat$fstatus == 1L) < 100L) return(list(coef = tibble(), curves = tibble()))

  sampled <- FALSE
  source_n <- nrow(source_dat)
  source_adverse <- sum(source_dat$fstatus == 1L)
  if (source_n > finegray_max_n) {
    set.seed(20260603 + match(org, target_organs) * 100 + match(model_name, names(model_specs)))
    sampled <- TRUE
    model_dat <- source_dat %>%
      group_by(fstatus) %>%
      slice_sample(prop = min(1, finegray_max_n / source_n)) %>%
      ungroup() %>%
      droplevels()
  } else {
    model_dat <- source_dat
  }
  model_dat <- model_dat %>%
    mutate(listing_center = collapse_top_levels(listing_center, n_top = 20L)) %>%
    droplevels()

  log_msg("Fine-Gray ", org, " | ", model_name, " n=", nrow(model_dat), " adverse=", sum(model_dat$fstatus == 1L))
  rhs <- c(exposure_terms, "age", "sex", "race", "index_year_centered", "organ_score", "listing_center")
  form <- as.formula(paste("~", paste(rhs, collapse = " + ")))
  pred_rows <- bind_rows(lapply(exposure_terms, function(focal) make_reference_rows(model_dat, exposure_terms, focal) %>% mutate(focal_term = focal)))
  combined <- bind_rows(
    model_dat %>% select(all_of(rhs)),
    pred_rows %>% select(all_of(rhs))
  )
  mm_all <- model.matrix(form, data = combined)
  mm_all <- mm_all[, colnames(mm_all) != "(Intercept)", drop = FALSE]
  mm_fit <- mm_all[seq_len(nrow(model_dat)), , drop = FALSE]
  mm_pred <- mm_all[(nrow(model_dat) + 1L):nrow(mm_all), , drop = FALSE]

  fit_error <- NULL
  elapsed <- system.time({
    fit <- tryCatch(
      crr(
        ftime = model_dat$followup_days,
        fstatus = model_dat$fstatus,
        cov1 = mm_fit,
        failcode = 1,
        cencode = 0
      ),
      error = function(e) {
        fit_error <<- conditionMessage(e)
        NULL
      }
    )
  })
  if (is.null(fit)) {
    coef_tbl <- tibble(
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      model = model_name,
      endpoint = "death_or_deterioration_delist_subdistribution",
      term = exposure_terms,
      exposure = recode(exposure_terms, !!!term_labels),
      source_n = source_n,
      source_adverse_events = source_adverse,
      n = nrow(model_dat),
      adverse_events = sum(model_dat$fstatus == 1L),
      sampled = sampled,
      finegray_max_n = finegray_max_n,
      subdistribution_hazard_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      elapsed_seconds = unname(elapsed[["elapsed"]]),
      center_adjustment = "top_20_listing_center_fixed_effect_dummies_plus_other",
      variance_estimator = "cmprsk_crr_model_based",
      error = fit_error
    )
    return(list(coef = coef_tbl, curves = tibble()))
  }

  coef <- fit$coef
  var <- diag(fit$var)
  keep <- intersect(exposure_terms, names(coef))
  coef_tbl <- tibble(
    organ = org,
    organ_label = recode(org, !!!organ_labels),
    model = model_name,
    endpoint = "death_or_deterioration_delist_subdistribution",
    term = keep,
    exposure = recode(keep, !!!term_labels),
    source_n = source_n,
    source_adverse_events = source_adverse,
    n = nrow(model_dat),
    adverse_events = sum(model_dat$fstatus == 1L),
    sampled = sampled,
    finegray_max_n = finegray_max_n,
    subdistribution_hazard_ratio = exp(coef[keep]),
    conf_low = exp(coef[keep] - 1.96 * sqrt(var[match(keep, names(coef))])),
    conf_high = exp(coef[keep] + 1.96 * sqrt(var[match(keep, names(coef))])),
    p_value = 2 * pnorm(abs(coef[keep] / sqrt(var[match(keep, names(coef))])), lower.tail = FALSE),
    elapsed_seconds = unname(elapsed[["elapsed"]]),
    center_adjustment = "top_20_listing_center_fixed_effect_dummies_plus_other",
    variance_estimator = "cmprsk_crr_model_based",
    error = NA_character_
  )

  pred <- predict(fit, cov1 = mm_pred)
  curve_tbl <- bind_rows(lapply(seq_len(nrow(mm_pred)), function(i) {
    tibble(
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      model = model_name,
      focal_term = pred_rows$focal_term[[i]],
      exposure = recode(pred_rows$focal_term[[i]], !!!term_labels),
      exposure_label = curve_terms$exposure_label[match(pred_rows$focal_term[[i]], curve_terms$term)],
      curve_group = pred_rows$curve_group[[i]],
      curve_exposure_value = pred_rows$curve_exposure_value[[i]],
      source_n = source_n,
      n = nrow(model_dat),
      sampled = sampled,
      finegray_max_n = finegray_max_n,
      time_days = pred[, 1],
      cif = pred[, i + 1L]
    )
  }))

  list(coef = coef_tbl, curves = curve_tbl)
}

coef_paths <- character()
curve_paths <- character()
for (org in target_organs) {
  dat <- base %>% filter(WL_ORG == org)
  for (model_name in names(model_specs)) {
    coef_path <- file.path(model_dir, paste0(tolower(org), "_", model_name, "_coefficients.csv"))
    curve_path <- file.path(model_dir, paste0(tolower(org), "_", model_name, "_curves.csv"))
    coef_paths <- c(coef_paths, coef_path)
    curve_paths <- c(curve_paths, curve_path)
    if (file.exists(coef_path) && file.exists(curve_path)) {
      log_msg("Skipping existing Fine-Gray model ", org, " | ", model_name)
      next
    }
    result <- fit_finegray_model(dat, org, model_name, model_specs[[model_name]])
    write_csv(result$coef, coef_path)
    write_csv(result$curves, curve_path)
    rm(result)
    invisible(gc())
  }
}

coef_results <- bind_rows(lapply(coef_paths[file.exists(coef_paths)], read_csv, show_col_types = FALSE))
curve_results <- bind_rows(lapply(curve_paths[file.exists(curve_paths)], read_csv, show_col_types = FALSE))
write_csv(coef_results, file.path(out_dir, "adjusted_finegray_all_exposure_coefficients.csv"))
write_csv(curve_results, file.path(out_dir, "adjusted_finegray_model_based_cif_curves.csv"))

plot_dat <- curve_results %>%
  filter(time_days <= 3650) %>%
  mutate(
    years = time_days / 365.25,
    curve_group = factor(curve_group, levels = c("Q1 reference", "Q2 reference", "Q4 reference"))
  )

plot_path <- file.path(fig_dir, "adjusted_finegray_model_based_cifs_all_exposures.png")
p <- ggplot(plot_dat, aes(x = years, y = cif, color = curve_group)) +
  geom_line(linewidth = 0.75) +
  facet_grid(exposure_label ~ organ_label, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_brewer(palette = "Dark2", na.translate = FALSE) +
  labs(
    title = "Adjusted Fine-Gray Model-Based CIFs by Exposure Level",
    subtitle = paste0(
      "Predicted death/deterioration delisting CIF at organ-specific Q1, median, and Q4 exposure values. ",
      "Models adjust for age, sex, race, listing year, organ score, and center fixed effects; capped at n=",
      finegray_max_n,
      " per organ/model when needed."
    ),
    x = "Years from waitlist start",
    y = "Predicted cumulative incidence",
    color = "Exposure value"
  ) +
  theme_cif()
ggsave(plot_path, p, width = 15, height = 11, dpi = 240, bg = "white")

write_csv(
  tibble(
    figure = "adjusted_finegray_model_based_cifs_all_exposures",
    path = plot_path
  ),
  file.path(fig_dir, "adjusted_finegray_model_based_cif_manifest.csv")
)

log_msg("Wrote adjusted Fine-Gray model-based CIF outputs to ", out_dir)
