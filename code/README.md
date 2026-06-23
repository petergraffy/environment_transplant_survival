# Code Directory

This directory contains the primary manuscript analysis scripts and shared
helpers. Scripts should be run from the repository root.

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
