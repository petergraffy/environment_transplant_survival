# Revised Waitlist Pollution Modeling Strategy

## Scientific Aim

Estimate whether ZIP/ZCTA-level air pollution adds independent prognostic
information for adverse waitlist outcomes after adjustment for existing
organ-specific clinical risk information, transplant center, calendar time, and
community-level confounding.

The revised analysis focuses on four organ systems:

- Heart (`HR`)
- Kidney (`KI`)
- Liver (`LI`)
- Lung (`LU`)

## Primary Outcome

Primary adverse waitlist outcome:

- `CAN_REM_CD == 8`: died
- `CAN_REM_CD == 13`: candidate condition deteriorated, too sick for transplant

Competing and descriptive outcome categories:

- Transplant/removal for transplant: `4`, `14`, `15`, `18`, `19`, `21`, `22`, `23`
- Improvement: `12`
- Medically unsuitable: `5`, retained for sensitivity analysis rather than mixed
  into the primary deterioration endpoint
- Other delisting/removal: other nonmissing codes

The SAF data dictionary defines `12` as condition improved and `13` as condition
deteriorated, too sick for transplant, so these should not be combined.

## Exposure Strategy

The revised model layer uses:

- Prior exposure: prior calendar year pollution at candidate ZIP. Because the
  pollution surface begins in 2005, the primary prior-exposure cohort begins in
  2006.
- Waitlist-duration exposure: lagged annual exposure updated over waitlist
  follow-up. Each interval in calendar year `Y` uses exposure year `Y - 1`, so
  exposure precedes the risk interval.
- Cumulative waitlist exposure: time-weighted mean of lagged annual exposures
  accrued before the start of each risk interval.

The annual time-updated approach is the first defensible implementation that
works across PM2.5, ozone, and annual NO2. A later refinement should use the
monthly PM2.5 and ozone surfaces for rolling 3-, 6-, and 12-month windows, while
retaining annual NO2.

## Community Covariates

`code/06_build_community_covariates.R` downloads ACS 5-year ZCTA covariates from
the Census API. The Census API currently requires `CENSUS_API_KEY` to be set in
the R session or shell before running the script.

- Median household income
- Poverty
- Bachelor's degree or higher
- Unemployment
- No household vehicle
- Nonwhite, Black, Hispanic composition

ACS 5-year ZCTA data are pulled for 2012-2023. Analysis years 2005-2011 are
backfilled to the 2012 ACS surface if needed.

Potential later additions:

- USDA RUCA ZIP code urbanicity
- CDC PLACES ZCTA health-behavior estimates
- ADI or SDI if tract/block-group data and a tract-to-ZCTA crosswalk are added

## Model Structure

The revised scripts fit organ-specific Cox models for adverse waitlist outcome.

Core adjustment:

- Age spline
- Sex, SRTR race, SRTR ethnicity
- BMI spline
- ABO
- Payor
- Education
- Community ACS covariates
- Calendar-year strata
- Transplant-center fixed effects via center strata
- Robust clustering by person

Organ-specific clinical adjustment:

- Heart: status, diagnosis, functional status, life support, ECMO, ventilator,
  IABP, inotropes, VAD/TAH fields, creatinine, albumin, diabetes.
- Kidney: dialysis, dialysis vintage, GFR, creatinine, PRA fields, diagnosis,
  functional status, diabetes.
- Liver: lab MELD, MELD type, bilirubin, creatinine, albumin, INR, sodium,
  ascites, encephalopathy, diagnosis, functional status, diabetes.
- Lung: status, diagnosis, functional status, life support, ECMO, ventilator,
  ICU, FVC, FEV1, PCO2, resting oxygen, six-minute walk, pulmonary pressures,
  cardiac output, corticosteroid dependence, diabetes.

## Scripts

- `code/06_build_community_covariates.R`: ACS ZCTA community covariates.
- `code/exploratory/07_organ_specific_adverse_waitlist_models.R`: revised heart/kidney/liver/lung
  organ-specific models with prior and waitlist-duration pollution exposure.

## Outputs

Revised outputs are written to:

- `output/organ_specific_adverse_waitlist/`

Key files:

- `revised_cohort_summary_by_organ.csv`
- `organ_specific_baseline_prior_pollution_cox_results.csv`
- `lagged_waitlist_exposure_interval_summary.csv`
- `organ_specific_prior_waitlist_pollution_timevarying_cox_results.csv`
