# September 2026 manuscript revision analyses

These scripts preserve the primary outputs and write to
`output/revision_2026_checklist`. Candidate-level intermediates are stored only
in its git-ignored `cache` directory. Run from the repository root with the
existing R environment and access to the Q1 2026 SAF and pollution inputs.

## Execution order

1. `01_audit_lung_sources.R` and `01b_audit_lung_offer_sources.R`: source audit.
2. `02_prepare_registration_sensitivities.R`: baseline cache, separate spells,
   and concurrent-registration intervals.
3. `03_baseline_diagnostics_eras.R`: exact primary reproduction, Schoenfeld
   diagnostics, follow-up bands, listing-era models and interaction tests.
4. `04_svi_validation.R`: observed 2022 ACS proxy versus official CDC SVI.
5. `05_aj_risk_tables.R`: all four AJ horizons with exact numbers at risk.
6. `06_recent_finegray.R`: 2020-2022 listings, 3-year cap, all eligible candidates.
7. `07_registration_models.R`: concurrency and separate-spell Cox sensitivities.
8. `08_timevarying_diagnostics.R`: exact primary reproduction and PH diagnostics
   for all 12 time-varying models; follow-up-band estimates.
9. `09_recent_cause_specific.R`: same-cohort, same-horizon Cox comparator.
10. `10_supplement_tables.R`: CSV, TSV and browser/Word-copyable HTML tables.
11. `11_validation_and_diagnostic_figures.R`: independent AJ-count checks and
    two-panel residual/coefficient-trend figures.
12. `12_visual_qa.py`: render all 34 PDFs and check paired PNG resolution and
    nonblank content; inspect the resulting contact sheets.
13. `13_thoracic_warning_audit.R`: independently audit thoracic/liver model
    warnings, recording the affected follow-up band rather than hiding them.
14. `14_baseline_warning_audit.R`: identify categorical coefficients responsible
    for warnings in baseline bands, policy eras and recent-cohort models.
15. `15_kidney_late_band_audit.R`: reproduce kidney late-band fits and capture
    their warnings, preserving the primary inclusive terminal-day convention.

Rerun step 10 after warning audits so the final tables include warning flags.

Example (PowerShell):

```powershell
& 'C:/Program Files/R/R-4.5.2/bin/Rscript.exe' code/tests/test_revision_design.R
& 'C:/Program Files/R/R-4.5.2/bin/Rscript.exe' code/revision_2026/07_registration_models.R
```

Run large time-varying diagnostics sequentially; do not run several copies.

`16_pre_covid_competing_risks.R` runs the updated 2015-2019 listing-cohort
Fine-Gray sensitivity and matched cause-specific models with 3-year follow-up.
It parameterizes scripts 06 and 09 via options and writes a separate
`finegray_2015_2019` directory, preserving the original 2020-2022 outputs.
After its validation passes, script 10 uses the updated cohort for the
supplemental sensitivity columns and exports a dedicated comparison table.
Pandemic-era follow-up is retained; only listings in 2020 onward are excluded.
The kidney models contain more than 20 million intervals. No random
subsampling or runtime timeout is used. Per-model diagnostic outputs and
per-organ/pollutant Fine-Gray checkpoints are written during the run.

The official SVI input is `data/reference/cdc_svi/SVI_2022_US_ZCTA.csv`, obtained
from https://svi.cdc.gov/Documents/Data/2022/csv/zcta/SVI_2022_US_ZCTA.csv.
Its recorded SHA256 is
`e2007404c6e63aca40fcfab64893599de8e3cd79887d61463b36b7ca60f2791a`.

## Interpretation safeguards

- No change is made to the primary adjustment sets, cohort or estimates.
- Concurrent count is updated using observed registrations, never future maximum
  concurrency. Gaps in the primary collapsed timeline receive a separate flag.
- Spells reset baseline covariates and exposure after a positive gap; overlapping
  and same-day touching registrations remain combined.
- Fine-Gray subdistribution HRs are not cause-specific HRs. Both center-specific
  censoring weights and center-stratified baseline hazards are used.
- Listing-era analyses are not causal estimates of policy effects.
- PH diagnostics can contradict the constant-coefficient assumption; report
  violations and time-band patterns rather than asserting automatic validity.
- Official lung LAS/CAS replacement is not implemented without a validated
  score-history source. See `docs/lung_score_source_audit.md`.
- SVI validation measures convergent validity, not equivalence. Median imputation
  of missing proxy components is retained and its frequency is reported.

Caches and Fine-Gray checkpoints belong to the frozen primary-input version
recorded in the output provenance. Use a new revision directory or explicitly
invalidate these caches after changing source data, cohort definitions or model
specifications; existing checkpoints are not automatically fingerprinted.

## Additional lung sensitivities

- `17_svi_validation_figure.R`: combined appendix SVI validation figure.
- `18_lung_score_reaudit.R`: fresh lung score/source schema audit.
- `19_inspect_las_inputs.R`, `20_lung_las_cas_sensitivity.R`, and
  `21_lung_las_cas_report.R`: exploratory reconstructed LAS/offer-level CAS
  sensitivity. The reconstructed LAS is default-dominated; do not treat it as
  validated clinical severity adjustment.
- `22_lung_composite_sensitivity.R` and `23_lung_composite_report.R`: fixed,
  study-defined clinical composite and pre/post-CAS era sensitivities, including
  single-/multipollutant and age/sex/race subgroup models. See
  `docs/lung_severity_composite_sensitivity.md` for weights and limitations.

Unit and output checks are in `code/tests/test_las_2021.R`,
`test_lung_score_timing.R`, `test_lung_composite.R`, and
`test_lung_composite_outputs.R`. All outputs remain separate from primary results.
