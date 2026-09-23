# Daily pollution inputs for the manuscript revision

The baseline and time-varying analysis scripts (78-81, 84, 86, 88, and 90)
use the published daily ZCTA releases through December 31, 2024:

- PM2.5: https://github.com/petergraffy/environment_transplant_survival/releases/tag/lghap-pm25-zcta-daily-v1
- O3: https://github.com/petergraffy/environment_transplant_survival/releases/tag/o3-zcta-daily-v1

Each release contains 20 annual Parquet files, 2005-2024, covering 33,300 ZCTAs.
PM2.5 is in ug/m3 and O3 is daily maximum 8-hour average concentration in ppb.
The local PM2.5 2022-2024 files are LGHAP estimates, not the separate local
AOD-derived experimental PM2.5 estimates. NO2 inputs are unchanged.

## Preparation and verification

Run from the repository root:

```powershell
& code/pollution_aggregation/verify_daily_analysis_releases.ps1
Rscript code/prepare_daily_pollution_inputs.R
Rscript code/prepare_rolling_prior_pollution.R
```

The verification script compares all local annual files with GitHub release
SHA-256 digests. The preparation script builds annual and monthly caches without
fitting clinical models, checks coverage and unique keys, and validates daily
averaging and leap-day inclusion against a sample ZCTA in 2024. Reports are
written to `output/exposure_input_validation/`.

The shared reader in `code/daily_pollution_inputs.R` requires every annual file,
processes one year at a time, and retains the existing minimum completeness
rules (300 days/year or 20 days/month). Cache names incorporate the source file
inventory, sizes, modification times, aggregation resolution, and column names.
Old unversioned caches are not reused. Changes to an input file trigger a new
cache. Run checksum verification after downloading or replacing release assets.

## Coverage and analysis implications

- Baseline Cox, AJ, subgroup, sensitivity, severity-association, and baseline
  map models now use listing-date-specific exposure from
  `code/rolling_prior_pollution.R`. PM2.5/O3 are the arithmetic mean of all 365
  daily values from listing date minus 365 through listing date minus 1.
  Listing-day exposure is excluded, and exactly 365 finite daily observations
  are required. Leap days count as ordinary days. Incomplete windows are missing;
  there is no shortened window or extrapolation. Current daily coverage supports
  listing dates January 1, 2006 through January 1, 2025, inclusive. Most 2025
  listings therefore lose PM2.5/O3 baseline eligibility compared with the former
  prior-calendar-year definition. Follow-up of eligible candidates can continue
  into 2025.
- NO2 uses the 12 complete calendar months before the listing month, weighted
  by each month's number of days. For example, a July 15, 2020 listing uses
  July 1, 2019 through June 30, 2020. All 12 months must have a finite exposure.
  The local monthly release begins January 2019 (not 2018); earlier months use
  the corresponding year's annual NO2 as an approximation, as requested.
  Annual values are never used to fill holes within the monthly-source era.
  Annual estimates cannot resolve within-year timing and may incorporate
  concentrations occurring after listing; this temporal limitation must be
  reported. `no2_prior_annual_months` records the number of approximated months.
  With these annual approximations, NO2 windows support listings from January
  2006 through the cohort's 2025 endpoint. Fully monthly windows begin January
  2020. NO2 windows exclude the entire listing month.
- Baseline outcome follow-up is not truncated at the exposure endpoint; exposure
  is fixed before listing and observed follow-up can continue into 2025.
- Time-varying PM2.5 follow-up now ends December 31, 2024 rather than December
  31, 2021. O3 already used December 31, 2024. Both use monthly means computed
  from the daily data, with no exposure carry-forward into 2025.
- Time-varying multipollutant models including PM2.5 now share coverage through
  December 31, 2024. NO2-only follow-up remains through December 31, 2025.
- The time-varying primary script now refits model results instead of silently
  reusing result files from the earlier exposure window.
- Concentration and survival-map input readers use the updated daily PM2.5
  data; the concentration-map reader also uses updated daily O3. Historical
  annual-input analyses in scripts 50, 63, and 72 remain unchanged.

This update prepares inputs only. Existing model estimates, figures, and tables
remain from previous runs until the corresponding analyses are explicitly rerun.
The survival-map estimand and other reviewer-requested methodological changes
are separate from this data-source update.

Rolling caches are keyed to the listing-date/ZCTA query set and source files,
separately from the annual caches. Window endpoints and valid-day/month counts
are retained for auditing. Preparation writes coverage by organ and listing
year to `output/exposure_input_validation/rolling_prior_exposure_coverage_by_organ_year.csv`.
The test script `code/tests/test_rolling_prior_pollution.R` checks exact window
boundaries, leap days, missing days, year partitions, duplicates, monthly
weighting, and annual-only NO2 fallback. National AJ quartiles will be
recomputed from the eligible candidates' rolling exposures when AJ is rerun.
