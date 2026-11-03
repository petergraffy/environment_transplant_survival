#!/usr/bin/env Rscript
# One-year RMST from the existing six Cox models; no clinical refitting required.
source("code/r_runtime.R"); ensure_user_library()
suppressPackageStartupMessages({library(data.table); library(survival); library(ggplot2)})
source("code/lung_survival_benefit_helpers.R")
setDTthreads(2)
dest <- Sys.getenv("LUNG_BENEFIT_OUTPUT", "output/lung_pollution_survival_benefit")
fitted <- readRDS(file.path(dest, "cache", "fitted.rds"))
horizon <- 365
exposure_curves <- list(); time_curves <- list()
for (p in c("pm25", "no2", "o3")) {
  message(format(Sys.time(), "%F %T"), " | One-year RMST / ", p)
  a <- fitted$models[[p]]$pretransplant; b <- fitted$models[[p]]$posttransplant
  r <- lung_supported_reference(fitted$reference[[p]], a, b, horizon)
  stopifnot(nrow(r) == nrow(fitted$reference[[p]]))
  h <- list(pre = lung_baseline(a), post = lung_baseline(b))
  exposure_curves[[p]] <- lung_predict_rmst(a, b, r, p, fitted$grids[[p]], horizon, h)
  time_curves[[p]] <- lung_predict_rmst(a, b, r, p, range(fitted$grids[[p]]),
    sort(unique(c(seq(0, 360, by = 30), horizon))), h)
}
curves <- rbindlist(exposure_curves); times <- rbindlist(time_curves)
fwrite(lung_validate_rmst(curves, times, fitted$point), file.path(dest, "rmst_validation.csv"))
curves[, confidence_status := "not_estimated"]
times[, confidence_status := "not_estimated"]
fwrite(curves, file.path(dest, "rmst_1yr_exposure_curves.csv"))
fwrite(times, file.path(dest, "rmst_1yr_time_curves.csv"))
contrast <- curves[order(concentration), .(low_exposure = first(concentration), high_exposure = last(concentration),
  rmst_benefit_low_days = first(rmst_benefit_days), rmst_benefit_high_days = last(rmst_benefit_days),
  benefit_difference_days = last(rmst_benefit_days) - first(rmst_benefit_days), recipients = first(recipients)), by = pollutant]
fwrite(contrast, file.path(dest, "rmst_1yr_exposure_contrasts.csv"))
labels <- c(pm25 = "PM2.5 (ug/m3)", no2 = "NO2 (ppb)", o3 = "O3 (ppb)")
for (d in list(curves, times)) d[, pollutant := factor(pollutant, levels = c("pm25", "no2", "o3"))]
has_ci <- "rmst_benefit_days_low" %in% names(curves) && any(curves$confidence_status == "estimated")
caption <- paste("Paired recipients; covariates fixed at the decision date. Observational model estimates.",
  if (has_ci) "Shading: pointwise 95% person-bootstrap confidence intervals." else "Point estimates; RMST confidence intervals pending.", sep = "\n")
benefit <- ggplot(curves, aes(concentration, rmst_benefit_days)) +
  (if (has_ci) geom_ribbon(aes(ymin = rmst_benefit_days_low, ymax = rmst_benefit_days_high),
    fill = "#226078", alpha = .18) else NULL) +
  geom_hline(yintercept = 0, color = "grey70") + geom_line(color = "#226078", linewidth = 1) +
  facet_wrap(~pollutant, scales = "free_x", labeller = as_labeller(labels)) +
  labs(title = "One-year lung transplant benefit in expected days alive",
    subtitle = "Difference in restricted mean survival time over 365 days",
    x = "Previous complete month's concentration", y = "Additional days alive with transplant", caption = caption) +
  theme_bw(base_size = 11)
component_metrics <- c("rmst_without_transplant_days", "rmst_after_transplant_days")
for (nm in c(paste0(component_metrics, "_low"), paste0(component_metrics, "_high")))
  if (!nm %in% names(curves)) curves[, (nm) := NA_real_]
components <- melt(curves, id.vars = c("pollutant", "concentration"),
  measure.vars = list(rmst_days = component_metrics, ci_low = paste0(component_metrics, "_low"),
    ci_high = paste0(component_metrics, "_high")), variable.name = "strategy")
components[, strategy := factor(strategy, levels = c(1, 2), labels = c("Without transplant", "After transplant"))]
component_plot <- ggplot(components, aes(concentration, rmst_days, color = strategy)) +
  (if (has_ci) geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = strategy), alpha = .15, colour = NA) else NULL) +
  geom_line(linewidth = 1) +
  facet_wrap(~pollutant, scales = "free_x", labeller = as_labeller(labels)) +
  scale_color_manual(values = c("#b86638", "#226078"), labels = c("Without transplant", "After transplant")) +
  scale_fill_manual(values = c("#b86638", "#226078"), guide = "none") +
  labs(title = "Expected days alive during the first year", x = "Previous complete month's concentration",
    y = "Restricted mean survival time (days)", color = NULL, caption = caption) +
  theme_bw(base_size = 11) + theme(legend.position = "right")
times[, exposure := factor(ifelse(concentration == min(concentration), "Lower exposure", "Higher exposure"),
  levels = c("Lower exposure", "Higher exposure")), by = pollutant]
time_plot <- ggplot(times, aes(horizon_days, rmst_benefit_days, color = exposure)) +
  (if (has_ci) geom_ribbon(aes(ymin = rmst_benefit_days_low, ymax = rmst_benefit_days_high, fill = exposure),
    alpha = .15, colour = NA) else NULL) +
  geom_hline(yintercept = 0, color = "grey70") + geom_line(linewidth = 1) +
  facet_wrap(~pollutant, labeller = as_labeller(labels)) +
  scale_color_manual(values = c("#b86638", "#226078")) +
  scale_fill_manual(values = c("#b86638", "#226078"), guide = "none") +
  scale_x_continuous(breaks = c(0, 90, 180, 270, 365)) +
  labs(title = "Accumulated transplant survival benefit through one year",
    subtitle = "Lower and higher exposure are the endpoints of each pollutant's evaluated range",
    x = "Days after the transplant decision", y = "Additional days alive with transplant", color = NULL, caption = caption) +
  theme_bw(base_size = 11) + theme(legend.position = "right")
for (name in c("benefit", "components", "time")) {
  fig <- switch(name, benefit = benefit, components = component_plot, time = time_plot)
  for (ext in c("png", "pdf")) ggsave(file.path(dest, paste0("rmst_1yr_", name, ".", ext)),
    fig, width = 12, height = 5, dpi = 170)
}
writeLines(c("RMST is the exact area under each fitted step-survival curve from 0 to 365 days, expressed in days alive.",
  "Before-transplant survival is conditional on surviving to the recipient's transplant decision date, with baseline hazard shifted by waiting time.",
  "The same recipient population, center strata, phase-specific ages, exposure scenarios and frozen LAS/covariates as the existing survival-benefit analysis are retained.",
  "RMST benefit = posttransplant RMST minus conditional no-transplant RMST. The difference can be negative when early operative risk outweighs early gains.",
  "Exposure curves report 365-day RMST across all 21 exposure scenarios. Time curves report cumulative RMST differences through 365 days at the lowest and highest evaluated exposures.",
  if (has_ci) "Pointwise 95% confidence intervals integrate RMST in joint person-bootstrap Cox refits; see rmst_bootstrap_methods.txt for details."
    else "These are observational model-based point estimates. Run code/94_lung_pollution_rmst_bootstrap.R for RMST confidence intervals."),
  file.path(dest, "rmst_methods.txt"))
print(contrast)
