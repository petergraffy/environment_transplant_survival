# Formal Waitlist Environment Analysis

## Primary Question

Estimate whether environmental exposures accrued during the waitlist period are associated with adverse waitlist outcomes among heart, kidney, liver, and lung candidates.

## Primary Waitlist-Period Air-Pollution Models

The primary air-pollution analysis is implemented in:

- `code/50_primary_waitlist_period_pollution_cox.R`

Outputs are written to:

- `output/primary_waitlist_period_pollution_cox/primary_waitlist_period_pollution_analysis_dataset.csv.gz`
- `output/primary_waitlist_period_pollution_cox/primary_waitlist_period_pollution_cox_results.csv`
- `output/primary_waitlist_period_pollution_cox/primary_waitlist_period_pollution_cohort_summary.csv`
- `output/primary_waitlist_period_pollution_cox/model_results/*.csv`

This analysis asks whether higher average air-pollution exposure experienced during the observed waitlist period is associated with adverse waitlist outcomes.

### Exposure Window

The exposure is a year/day-weighted mean over the candidate's observed waitlist period:

- `PM2.5`: annualized from monthly ZCTA PM2.5 and followed through 2023
- `O3`: annualized from monthly ZCTA ozone and followed through 2023
- `NO2`: annual ZCTA NO2 and followed through 2025

Candidates still on the waitlist after the pollutant-specific exposure-data end year are censored at the end of that exposure period for that pollutant model.

### Primary Model Form

```r
Surv(followup_days_pollutant, pollutant_specific_adverse_event) ~
  waitlist_period_pollution +
  age + sex + race + listing_year + organ_score +
  strata(listing_center)
```

This is a cause-specific Cox model with model-based standard errors. The adverse event is death or delisting due to deterioration; transplant, improvement delisting, other removal, and administrative end of exposure follow-up are treated as censoring events.

## Baseline Annual Air-Pollution Study

The baseline listing-year air-pollution models are implemented in:

- `code/49_baseline_air_pollution_waitlist_cox.R`

Outputs are written to:

- `output/baseline_air_pollution_waitlist_cox/baseline_air_pollution_waitlist_analysis_dataset.csv.gz`
- `output/baseline_air_pollution_waitlist_cox/baseline_air_pollution_cox_results.csv`
- `output/baseline_air_pollution_waitlist_cox/baseline_air_pollution_cohort_summary.csv`
- `output/baseline_air_pollution_waitlist_cox/model_results/*.csv`

### Cohort

- SAF release: Q1 2026 SRTR SAF, resolved from `C:/Users/Peter Graffy/Box/Q1 2026 SAF` or `SAF_Q1_2026_DIR`
- Organs: heart (`HR`), kidney (`KI`), liver (`LI`), lung (`LU`)
- Listing years: 2006-2023 by default, matching the period with full PM2.5 and O3 exposure data
- Candidate sources: active/relisted candidates (`CAN_SOURCE %in% c("A", "R")`)
- Follow-up: activation/listing date to death, waitlist removal, last follow-up, or administrative censoring at 2023-12-31

### Outcome

The endpoint is death or delisting due to deterioration:

- adverse event: `CAN_REM_CD %in% c(8, 13)`
- successful transplant and delisting due to improvement are censored in the cause-specific Cox model
- other non-adverse removals and administrative end of follow-up are also treated as censoring events

### Exposures

Exposures are annual ZCTA values at baseline, defined by the candidate residential ZCTA and listing year.

- `NO2`: annual ZCTA NO2 for the listing year
- `PM2.5`: annual mean of monthly ZCTA PM2.5 for the listing year
- `O3`: annual mean of monthly ZCTA ozone for the listing year

The script now fits two exposure windows:

- `listing_year`: the listing-year value only
- `listing_prior_mean`: the mean of the listing-year and prior-year ZCTA values

No PM2.5 or O3 values are carried forward beyond 2023 in this baseline analysis.

Models are fit independently for each exposure and organ.

Scaling:

- `no2_10unit`: per 10-unit increment
- `pm25_5ug`: per 5 ug/m3
- `o3_10ppb`: per 10 ppb

### Model Form

```r
Surv(followup_days, adverse_event) ~
  exposure +
  age + sex + race + listing_year + organ_score +
  strata(listing_center)
```

This is a cause-specific Cox model with listing-center stratified baseline hazards and model-based standard errors.

## Primary Static Exposure Models

The completed primary static-exposure models are implemented in:

- `code/31_formal_waitlist_environment_primary_cox.R`

Outputs are written to:

- `output/formal_waitlist_environment/formal_waitlist_environment_analysis_dataset.csv.gz`
- `output/formal_waitlist_environment/formal_primary_cox_exposure_results.csv`
- `output/formal_waitlist_environment/formal_primary_cox_cohort_summary.csv`
- `output/formal_waitlist_environment/primary_cox_model_results/*.csv`

### Cohort

- Organs: heart (`HR`), kidney (`KI`), liver (`LI`), lung (`LU`)
- Listing years: 2006-2023
- Candidate sources: active/relisted candidates (`CAN_SOURCE %in% c("A", "R")`)
- Follow-up: activation/listing date to observed removal/death/end follow-up, administratively censored at 2023-12-31

### Outcome

The primary endpoint is death or delisting for deterioration:

- adverse event: `CAN_REM_CD %in% c(8, 13)`
- transplant, improvement delisting, other removal, or administrative censoring are treated as non-events in the cause-specific Cox model

### Exposures

All exposures are averaged from waitlist start through event/censor date.

- `tmax`: exact mean of daily gridMET maximum temperature over observed waitlist days
- `rmax`: exact mean of daily gridMET maximum relative humidity over observed waitlist days
- `PM2.5`: waitlist-day-weighted mean of monthly ZCTA PM2.5
- `O3`: waitlist-day-weighted mean of monthly ZCTA ozone
- `NO2`: waitlist-day-weighted mean of annual ZCTA NO2

Scaling:

- `tmax_waitlist_5c`: per 5 C
- `rmax_waitlist_10pct`: per 10 percentage points
- `pm25_waitlist_5ug`: per 5 ug/m3
- `o3_waitlist_10ppb`: per 10 ppb
- `no2_waitlist_10unit`: per 10-unit increment

### Model Sequence

Models are fit separately by organ. Exposure models are independent except heat/humidity:

- heat/humidity model: `tmax + rmax`
- PM2.5 model: `pm25`
- ozone model: `o3`
- NO2 model: `no2`

Primary model form:

```r
Surv(followup_days, adverse_event) ~
  exposure_terms +
  age + sex + race + index_year_centered + organ_score +
  listing_center +
  cluster(PERS_ID)
```

This is a cause-specific Cox model with listing-center fixed-effect indicators and robust clustering by patient.

### Baseline Organ Score

- Heart: baseline US-CRS proxy from the first heart status-justification record, with candidate-table fallback and median imputation for missing heart score labs
- Kidney: baseline dialysis + age proxy
- Liver: initial SRTR MELD/PELD, with last SRTR MELD/PELD used only if initial is missing
- Lung: baseline LAS/CAS component proxy from available candidate-table clinical fields

## Fine-Gray Note

We attempted full Fine-Gray subdistribution models with listing-center fixed effects. Both `survival::finegray()` plus weighted Cox and direct `cmprsk::crr()` were too slow at this sample size with center fixed-effect indicators. Because `coxph()` directly supports robust patient clustering and center fixed effects at this scale, the current formal primary model is cause-specific Cox.

A reduced Fine-Gray sensitivity can still be added, but it will likely need one of these compromises:

- center stratification rather than explicit fixed-effect coefficients
- no center fixed effects
- organ-specific sampling for computational benchmarking
- high-performance implementation outside base R/cmprsk

## Time-Varying Organ Score Feasibility

The SAF status-history files do not provide the same longitudinal score components for every organ.

- Heart: updated US-CRS components are available in `statjust_hr1a.sas7bdat` by `CANHX_CHG_DT`; these can be joined to follow-up intervals as most recent prior status-justification values.
- Liver: updated MELD information is available in `stathist_liin.sas7bdat`, including `CANHX_SRTR_LAB_MELD`, `CANHX_OPTN_LAB_MELD`, and related MELD laboratory fields.
- Kidney: `stathist_kipa.sas7bdat` contains status intervals and CPRA, but not the full dialysis + age score components beyond candidate-level dialysis dates. The feasible dynamic proxy is updated age and dialysis duration over time.
- Lung: `stathist_thor.sas7bdat` contains status intervals but not longitudinal LAS/CAS component fields. The feasible model can update waitlist status over time, but the LAS/CAS component proxy remains baseline unless another longitudinal lung-status source is identified.

The next time-varying Cox model should therefore be organ-specific rather than pretending all four scores are equally updateable.

## Time-Varying Cox Models

The time-varying interval build and model fitters are:

- `code/32_timevarying_cox_and_finegray_sensitivity.R`
- `code/33_fit_timevarying_models_from_intervals.R`

Outputs are written to:

- `output/formal_waitlist_environment_timevarying/timevarying_interval_analysis_dataset.csv.gz`
- `output/formal_waitlist_environment_timevarying/timevarying_cox_exposure_results.csv`
- `output/formal_waitlist_environment_timevarying/timevarying_cox_model_results/*.csv`

The interval builder splits each candidate's observed waitlist follow-up by SAF status-history intervals. Environmental exposures are averaged within each interval using month/day weights:

- `tmax` and `rmax`: ZCTA-month gridMET means weighted by interval days in each month
- `PM2.5` and `O3`: ZCTA-month pollution means weighted by interval days
- `NO2`: annual ZCTA means weighted by interval days

Time-varying organ score handling:

- Heart: most recent prior heart status-justification US-CRS proxy where available; baseline score otherwise
- Kidney: updated age and dialysis-duration proxy over interval time
- Liver: interval SRTR/OPTN MELD from liver status history where available; baseline score otherwise
- Lung: baseline LAS/CAS component proxy because longitudinal lung score components are not present in status history

Model form:

```r
Surv(tstart, tstop, tv_adverse_event) ~
  exposure_terms +
  age_interval + sex + race + index_year_centered + organ_score_tv +
  strata(listing_center) +
  cluster(PERS_ID)
```

Center is handled by center-stratified baseline hazards rather than center dummy coefficients for tractability with millions of intervals.

## Reduced Fine-Gray Heat/Humidity Sensitivity

The Fine-Gray heat/humidity sensitivity is implemented in:

- `code/34_reduced_finegray_heat_humidity_month_weighted.R`

Output:

- `output/formal_waitlist_environment_timevarying/reduced_finegray_heat_humidity_results.csv`

This uses the month/day-weighted waitlist-period `tmax` and `rmax` values from `output/waitlist_environment_exposure/waitlist_environment_exposures.csv.gz`. To make `cmprsk::crr()` computationally practical, it uses a stratified sample capped by `FINEGRAY_HEAT_MAX_N`; the completed run used 10,000 records per organ. It is therefore a sensitivity/diagnostic analysis rather than a definitive full-cohort Fine-Gray model.
