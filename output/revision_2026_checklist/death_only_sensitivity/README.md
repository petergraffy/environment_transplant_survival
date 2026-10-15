# Death-only sensitivity analysis

Values are adjusted cause-specific hazard ratios (95% CIs). Models use the identical cohort, exposure, covariates, and follow-up as the corresponding primary model; only the endpoint changes.
Baseline: exposure in the 365 days before listing for PM2.5 and O3, or 12 complete calendar months before listing for NO2 (annual approximations when monthly data are unavailable); adjustment for age, sex, race, ACS-derived ZCTA vulnerability proxy, and listing-center strata.
Time-varying: monthly pollution intervals and the existing primary severity covariates, in addition to age, sex, race, vulnerability proxy, and center strata. Kidney: no-dialysis indicator, updated dialysis duration, and diabetes. Heart/liver: existing updated score; lung: existing primary component proxy. This analysis does not substitute the separate lung-composite sensitivity score. Neither model includes listing year.
Death-only events are final waitlist removals coded 8 (died). Final deterioration removals (code 13), transplantation, improvement, and other removals are censored at their recorded exit dates; administrative and pollutant-coverage censoring are unchanged. The primary registration-collapsing rules are retained. Death after final removal is not counted. Transplant-procedure deaths (codes 21/23) retain their primary transplant classification.
These are cause-specific hazard estimates, not cumulative probabilities. Censoring deterioration removals does not estimate mortality if those candidates had remained on the waitlist. Model-based standard errors and Efron ties were used; no person-level robust clustering.
[a] Four lung PM2.5/O3 death-only fits had a separation warning for the missing-race coefficient, not the pollutant coefficient. No deaths occurred in that race category. Diagnostic refits excluding only that category changed the pollutant HR by less than 0.0001 and eliminated the warnings. Full-cohort estimates are retained here; diagnostic refits are not substituted. See lung_warning_audit.csv.
Abbreviations: ACS, American Community Survey; CI, confidence interval; HR, hazard ratio; NO2, nitrogen dioxide; O3, ozone; PM2.5, particulate matter with diameter of 2.5 micrometers or less; ZCTA, ZIP Code Tabulation Area.

## Validation
All 24 composite refits matched the saved primary HRs within 0.00001, with identical candidate and composite-event counts. Outcome components reconcile exactly to the full primary cohort.
Death-only point estimates had the same direction as the corresponding composite model in 18 of 24 comparisons.

## Reproduction
Run code/revision_2026/25_death_only_sensitivity.R, code/revision_2026/27_death_only_warning_audit.R, then code/revision_2026/26_death_only_report.R. Existing model-pair outputs are checkpoints; remove them before rerunning with changed inputs. Run code/figures_tables/69_make_table1_baseline_characteristics.R to regenerate the full Table 1.

## Files
death_only_sensitivity.html: formatted supplemental comparison table.
death_only_comparison.csv: numeric comparison and formatted HRs.
full_model_results.csv: all log-HR coefficients, SEs, CIs, and fit warnings.
outcome_components_by_organ.csv / removal_code_audit.csv: aggregate outcome coding audit.
Individual-level outcome linkage is kept only in the ignored private cache.
