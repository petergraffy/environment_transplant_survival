# Environment Transplant Survival

This repository contains the analysis code and generated display assets for a
study of air pollution exposure and adverse outcomes among heart, kidney, liver,
and lung transplant waitlist candidates.

## Study Overview

The primary analysis estimates associations between waitlist-period air
pollution exposure and death or delisting due to clinical deterioration using
organ-specific cause-specific Cox proportional hazards models. Exposures are
day-weighted mean pollutant concentrations during observed waitlist follow-up:

- PM2.5, scaled per 5 ug/m3
- NO2, scaled per 10 ppb
- O3, scaled per 10 ppb

Models are fit separately by organ and pollutant. Heart, liver, and lung models
adjust for age, sex, race, listing year, organ-specific baseline severity, and
listing-center strata. Kidney models adjust for age, sex, race, listing year,
no dialysis time, dialysis duration, diabetes, and listing-center strata.

## Repository Layout

- `code/`: primary analysis scripts and shared helpers.
- `code/figures_tables/`: scripts that create manuscript figures,
  supplemental figures, and tables from analysis outputs.
- `code/pollution_aggregation/`: upstream scripts that aggregate, QA, package,
  and upload ZCTA-level air pollution exposure surfaces.
- `code/exploratory/`: earlier exploratory, diagnostic, time-varying,
  competing-risk, gridMET, and data-structure scripts not required for the
  primary manuscript analysis.
- `data/release/`: derived exposure release assets used by the analysis.
- `docs/`: notes describing exposure-surface construction and release assets.
- `output/`: generated model outputs, figures, and tables.

## Primary Analysis Scripts

Run scripts from the repository root so relative paths resolve correctly.

1. `code/50_primary_waitlist_period_pollution_cox.R`
   - Builds the deduplicated waitlist cohort.
   - Collapses same-person, same-organ registrations to first-listing baseline
     exposure/covariates and final observed waitlist outcome.
   - Computes waitlist-period pollutant exposures.
   - Fits the primary cause-specific Cox models.

2. `code/63_primary_waitlist_pollution_subgroup_cox.R`
   - Fits stratified cause-specific Cox models by age group, sex, and race.

3. `code/72_primary_waitlist_pollution_cox_sensitivity_models.R`
   - Fits sensitivity models, including models without organ score, without
     listing-center strata, with an ACS-derived ZCTA social vulnerability proxy,
     prior-year exposure models, and multipollutant models.

4. `code/figures_tables/`
   - Generates manuscript and supplemental figures/tables from the outputs of
     the model scripts.

## Data Requirements

The analysis uses SRTR Standard Analysis Files (SAF), local ZCTA-level pollution
exposure assets, and ACS-derived ZCTA community covariates. SRTR SAF files are
not included in this repository. Local SAF paths are configured in
`code/saf_paths.R`; helper runtime setup is in `code/r_runtime.R`.

## Reproducibility Notes

- Most scripts are designed to be run with `Rscript` from the repository root.
- Output files are written under `output/`.
- Large external raw data files and restricted SRTR source data are expected to
  exist locally and are not redistributed here.
- The committed outputs reflect the analysis state used for manuscript drafting.
