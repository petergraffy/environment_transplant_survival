# Code Directory

This directory contains the primary manuscript analysis scripts and shared
helpers. Scripts should be run from the repository root.

For the revised baseline and time-varying analyses, daily PM2.5 and O3 now use
the published 2005-2024 releases. Run `prepare_daily_pollution_inputs.R` to
prepare and validate exposures without fitting models. Source details,
coverage implications, and checksum verification are documented in
[`docs/daily_pollution_analysis_inputs.md`](../docs/daily_pollution_analysis_inputs.md).

Primary scripts:

- `50_primary_waitlist_period_pollution_cox.R`: builds the analysis cohort,
  computes waitlist-period pollution exposures, and fits primary cause-specific
  Cox models.
- `63_primary_waitlist_pollution_subgroup_cox.R`: fits stratified Cox models.
- `72_primary_waitlist_pollution_cox_sensitivity_models.R`: fits sensitivity
  models, including multipollutant models.
- `06_build_community_covariates.R`: builds ACS-derived ZCTA community
  covariates used for the social vulnerability proxy sensitivity analysis.
- `saf_paths.R` and `r_runtime.R`: local path and runtime helpers.

Subfolders:

- `figures_tables/`: figure and table creation scripts.
- `pollution_aggregation/`: upstream pollution exposure aggregation scripts.
- `exploratory/`: exploratory, diagnostic, and legacy analysis scripts.
