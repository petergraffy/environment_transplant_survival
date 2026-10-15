source("code/revision_2026/common.R")
library(gt)
dest <- file.path(revision_dir, "death_only_sensitivity")
r <- fread(file.path(dest, "pollutant_results.csv"))
r[, warnings := fifelse(is.na(warnings), "", as.character(warnings))]
stopifnot(nrow(r) == 48, all(is.finite(r$hazard_ratio)), all(r$HR_low < r$HR_high))
r[, effect := sprintf("%.2f (%.2f-%.2f)", hazard_ratio, HR_low, HR_high)]
wide <- merge(
  r[endpoint == "Composite", .(model, organ, pollutant, N = candidate_organ_episodes,
    composite_events = events, composite_HR = effect, composite_logHR = estimate)],
  r[endpoint == "Death only", .(model, organ, pollutant, death_N = candidate_organ_episodes,
    death_events = events, death_HR = effect, death_p = p.value, death_logHR = estimate, warnings)],
  by = c("model", "organ", "pollutant"))
stopifnot(all(wide$N == wide$death_N), all(wide$death_events <= wide$composite_events))
wide[, `:=`(same_direction = sign(composite_logHR) == sign(death_logHR),
  Organ = factor(organ_names[organ], levels = unname(organ_names)),
  Pollutant = factor(pollutant, levels = c("pm25", "no2", "o3"),
    labels = c("PM2.5 (per 5 ug/m3)", "NO2 (per 10 ppb)", "O3 (per 10 ppb)")))]
wide[nzchar(warnings), death_HR := paste0(death_HR, " [a]")]
setorder(wide, model, Organ, Pollutant)
fwrite(wide, file.path(dest, "death_only_comparison.csv"))
notes <- c(
  "Values are adjusted cause-specific hazard ratios (95% CIs). Models use the identical cohort, exposure, covariates, and follow-up as the corresponding primary model; only the endpoint changes.",
  "Baseline: exposure in the 365 days before listing for PM2.5 and O3, or 12 complete calendar months before listing for NO2 (annual approximations when monthly data are unavailable); adjustment for age, sex, race, ACS-derived ZCTA vulnerability proxy, and listing-center strata.",
  "Time-varying: monthly pollution intervals and the existing primary severity covariates, in addition to age, sex, race, vulnerability proxy, and center strata. Kidney: no-dialysis indicator, updated dialysis duration, and diabetes. Heart/liver: existing updated score; lung: existing primary component proxy. This analysis does not substitute the separate lung-composite sensitivity score. Neither model includes listing year.",
  "Death-only events are final waitlist removals coded 8 (died). Final deterioration removals (code 13), transplantation, improvement, and other removals are censored at their recorded exit dates; administrative and pollutant-coverage censoring are unchanged. The primary registration-collapsing rules are retained. Death after final removal is not counted. Transplant-procedure deaths (codes 21/23) retain their primary transplant classification.",
  "These are cause-specific hazard estimates, not cumulative probabilities. Censoring deterioration removals does not estimate mortality if those candidates had remained on the waitlist. Model-based standard errors and Efron ties were used; no person-level robust clustering.",
  "[a] Four lung PM2.5/O3 death-only fits had a separation warning for the missing-race coefficient, not the pollutant coefficient. No deaths occurred in that race category. Diagnostic refits excluding only that category changed the pollutant HR by less than 0.0001 and eliminated the warnings. Full-cohort estimates are retained here; diagnostic refits are not substituted. See lung_warning_audit.csv.",
  "Abbreviations: ACS, American Community Survey; CI, confidence interval; HR, hazard ratio; NO2, nitrogen dioxide; O3, ozone; PM2.5, particulate matter with diameter of 2.5 micrometers or less; ZCTA, ZIP Code Tabulation Area."
)
display <- as.data.frame(wide[, .(model, Organ, Pollutant, N, composite_events, composite_HR, death_events, death_HR)])
tab <- gt(display, groupname_col = "model") |>
  tab_header(title = "Death-Only Sensitivity Analysis") |>
  text_transform(locations = cells_body(columns = Pollutant), fn = function(x) {
    x <- gsub("PM2.5", "PM<sub>2.5</sub>", x, fixed = TRUE)
    x <- gsub("NO2", "NO<sub>2</sub>", x, fixed = TRUE)
    x <- gsub("O3", "O<sub>3</sub>", x, fixed = TRUE)
    lapply(gsub("ug/m3", "&#181;g/m<sup>3</sup>", x, fixed = TRUE), html)
  }) |>
  cols_label(N = "Candidates", composite_events = "Death/deterioration events",
    composite_HR = "Primary HR (95% CI)", death_events = "Deaths",
    death_HR = "Death-only HR (95% CI)") |>
  fmt_integer(columns = c(N, composite_events, death_events), use_seps = TRUE) |>
  tab_options(table.font.names = "Arial", table.font.size = px(14),
    column_labels.font.weight = "bold", row_group.font.weight = "bold", table.width = pct(100))
for (note in notes) tab <- tab |> tab_source_note(note)
gtsave(tab, file.path(dest, "death_only_sensitivity.html"))

counts <- fread(file.path(dest, "outcome_components_by_organ.csv"))
labels <- c(death = "Death on the waitlist", deterioration = "Removal for deterioration",
  transplant = "Transplantation", improvement = "Removal for improvement",
  other_exit = "Other waitlist exit", administrative_censor = "Administrative censoring")
long <- melt(counts, id.vars = c("WL_ORG", "N"), measure.vars = names(labels), variable.name = "component", value.name = "n")
long[, `:=`(Organ = organ_names[WL_ORG], Outcome = factor(labels[as.character(component)], levels = unname(labels)),
  value = sprintf("%s (%.1f%%)", format(n, big.mark = ",", trim = TRUE), 100 * n / N))]
ct <- dcast(long, Outcome ~ Organ, value.var = "value")
setcolorder(ct, c("Outcome", unname(organ_names)))
fwrite(ct, file.path(dest, "table1_outcome_components.csv"))
ctab <- gt(as.data.frame(ct)) |>
  tab_header(title = "Components of Waitlist Outcomes by Organ") |>
  tab_source_note("Values are n (% of the full organ-specific cohort). Outcomes are mutually exclusive, before pollutant-specific censoring. Death plus deterioration is the primary adverse endpoint; transplant plus improvement is the competing endpoint.") |>
  tab_source_note("Transplantation includes deaths during the transplant procedure, preserving the primary outcome definition: heart 85, kidney 10, liver 542, and lung 33.") |>
  tab_options(table.font.names = "Arial", table.font.size = px(14), column_labels.font.weight = "bold")
gtsave(ctab, file.path(dest, "table1_outcome_components.html"))
writeLines(c("# Death-only sensitivity analysis", "", notes,
  "", "## Validation", "All 24 composite refits matched the saved primary HRs within 0.00001, with identical candidate and composite-event counts. Outcome components reconcile exactly to the full primary cohort.",
  paste0("Death-only point estimates had the same direction as the corresponding composite model in ", sum(wide$same_direction), " of 24 comparisons."),
  "", "## Reproduction", "Run code/revision_2026/25_death_only_sensitivity.R, code/revision_2026/27_death_only_warning_audit.R, then code/revision_2026/26_death_only_report.R. Existing model-pair outputs are checkpoints; remove them before rerunning with changed inputs. Run code/figures_tables/69_make_table1_baseline_characteristics.R to regenerate the full Table 1.",
  "", "## Files", "death_only_sensitivity.html: formatted supplemental comparison table.",
  "death_only_comparison.csv: numeric comparison and formatted HRs.",
  "full_model_results.csv: all log-HR coefficients, SEs, CIs, and fit warnings.",
  "outcome_components_by_organ.csv / removal_code_audit.csv: aggregate outcome coding audit.",
  "Individual-level outcome linkage is kept only in the ignored private cache."), file.path(dest, "README.md"))
print(wide[, .(model, Organ, Pollutant, composite_HR, death_HR, warnings)])
