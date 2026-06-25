# Score Benchmark and Pollution Models

The benchmark-analysis layer asks whether air pollution improves prediction of
adverse waitlist outcome beyond the organ-specific score or score proxy.

Primary outcome:

- Death on the waitlist (`CAN_REM_CD == 8`)
- Delisted because condition deteriorated, too sick for transplant (`CAN_REM_CD == 13`)

Organ benchmark definitions:

- Heart: US-CRS proxy from heart status-justification variables. Inputs include
  short-term MCS, durable LVAD, bilirubin, race-neutral CKD-EPI 2021 eGFR,
  albumin, sodium, and BNP. The current implementation uses the corrected
  US-CRS formula:
  `1.02*short_MCS + 0.55*log(bilirubin+1) - 0.01*eGFR +
  0.40*log(BNP+1)*BNP_indicator + 0.20*log(BNP+1)*NT_proBNP_indicator -
  0.63*albumin - 0.07*sodium - 1.12*durable_LVAD`. BNP type is taken from
  `RiskStratDataHR` when available; statjust-only BNP values without type are
  coded as unknown type and assigned the BNP coefficient. Missing numeric
  components are median-imputed among heart candidates.
- Liver: initial SRTR lab MELD, falling back to last SRTR lab MELD.
- Lung: LAS/CAS component proxy from available thoracic candidate physiology and
  support variables because no direct LAS/CAS field was found in the SAF headers
  inspected. Replace this proxy with direct LAS/CAS values if they are found in
  another SAF supplement or match-run file.
- Kidney: dialysis plus age benchmark, using age, dialysis status, and dialysis
  vintage.

Models compared:

- Score only
- Score plus center fixed effects, calendar year, and ACS ZCTA community
  covariates
- Score/context plus prior-year pollution
- Score/context plus prior-year and current-year pollution

Script:

- `code/exploratory/08_score_benchmark_pollution_models.R`

Outputs:

- `output/score_benchmark_pollution/score_benchmark_cohort_summary.csv`
- `output/score_benchmark_pollution/score_benchmark_model_performance.csv`
- `output/score_benchmark_pollution/score_benchmark_pollution_coefficients.csv`
- `output/score_benchmark_pollution/score_benchmark_incremental_performance.csv`
- `output/score_benchmark_pollution/heart_uscrs_42day_model_performance.csv`
- `output/score_benchmark_pollution/heart_uscrs_42day_pollution_coefficients.csv`
- `output/score_benchmark_pollution/heart_uscrs_42day_incremental_performance.csv`

Current caveats:

- These are setup/benchmark outputs, not final inferential models.
- Heart US-CRS should be evaluated on its intended short horizon; the script now
  includes a 42-day adverse waitlist endpoint for heart in addition to full
  follow-up outputs.
- The current-year pollution terms are annual index-year exposures, not the full
  lagged cumulative waitlist exposure intervals from
  `code/exploratory/07_organ_specific_adverse_waitlist_models.R`.
- Heart and lung benchmarks need score-specific refinement before manuscript
  use: direct LAS/CAS recovery or a validated proxy for lung.
