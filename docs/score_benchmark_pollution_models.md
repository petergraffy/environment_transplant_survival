# Score Benchmark and Pollution Models

The benchmark-analysis layer asks whether air pollution improves prediction of
adverse waitlist outcome beyond the organ-specific score or score proxy.

Primary outcome:

- Death on the waitlist (`CAN_REM_CD == 8`)
- Delisted because condition deteriorated, too sick for transplant (`CAN_REM_CD == 13`)

Organ benchmark definitions:

- Heart: US-CRS proxy from heart status-justification variables. Inputs include
  short-term MCS, durable LVAD, bilirubin, race-neutral CKD-EPI 2021 eGFR,
  albumin, sodium, and BNP. The current implementation uses the US-CRS log-odds
  formula:
  `-0.656*albumin + 0.617*log(bilirubin+1) - 0.012*eGFR - 0.077*sodium -
  0.377*LVAD + 1.092*short_term_MCS + 0.433*log(BNP+1)`. The SAF field found so
  far is `CANHX_LAB_BNP`, without a BNP vs NT-proBNP type indicator, so this
  implementation currently treats it as regular BNP. Missing numeric components
  are median-imputed among heart candidates.
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
  use: BNP vs NT-proBNP handling for heart if a type field is found, and direct
  LAS/CAS recovery or a validated proxy for lung.
