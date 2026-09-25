# Lung score source audit

Audit date: September 23, 2026. Source: Q1 2026 SRTR SAF.

## Finding

The current score-adjusted models do **not** use separate official LAS and CAS
scores. `code/50_primary_waitlist_period_pollution_cox.R` constructs
`lung_las_cas_component_score`, a weighted clinical-component proxy.
`code/79_timevarying_pollution_cox_svi.R` retains this baseline lung proxy while
updating heart and liver severity inputs. Table 1 correctly calls the lung
measure a proxy linear predictor; it should not be described as official LAS,
CAS, or longitudinally updated lung urgency.

## Source inspection

All non-PTR SAS file headers in the available SAF were inspected for lung
allocation, composite allocation, urgency, LAS and CAS fields. No candidate-level
official LAS/CAS score-history source was identified. Thoracic status codes were
also inspected; they are status categories, not continuous official scores.

The 2018-2025 lung PTR offer-file headers contain `PTR_TOT_SCORE` and
`MATCH_SUBMIT_DT`. The accompanying `PtrLUFileFormat.pdf` defines these as
"Total score for the candidate on the match" and the match submission date/time.
The dictionary does not establish this field as the desired longitudinal
medical-urgency component. Observations are conditional on appearing in a
match run, not on every listing or clinical-score update. A first later offer
score cannot safely be substituted for a baseline score.

Machine-readable audit outputs are in `output/revision_2026_checklist/`:
`saf_variable_schema.csv`, `lung_score_source_candidates.csv`,
`lung_offer_variable_schema.csv`, and `lung_history_status_codes.csv`.

## September 24 value-level re-audit

Fresh headers were read for every non-PTR SAS table and all eight annual lung
PTR files in `ptr_lu_2018_2025`. The local data dictionary was also searched
for LAS, CAS, allocation scores and medical urgency. No explicit candidate-level
LAS, CAS medical-urgency or WLAUC history field was found.

Unlike the earlier header-only PTR inspection, the re-audit read score values:

- `PTR_TOT_SCORE` is missing for every record in 2018-2022 and all 346,544
  pre-March-9 records in the 2023 file. It cannot supply historical LAS here.
- In the CAS portion of 2023, 2,899 of 1,364,751 records lack the score.
  All 1,713,739 records in 2024 and 1,659,476 in 2025 have populated scores.
- In 2024, 251,509 of 285,577 registration-days with scores have multiple
  different scores. This is consistent with donor-specific allocation scoring;
  it cannot be assumed to represent a change in clinical severity.

These are counts of match records/registration-days, not cohort candidate counts.
Outputs and reproducible code: `output/revision_2026_checklist/lung_score_reaudit/`
and `code/revision_2026/18_lung_score_reaudit.R`.

Official documentation confirms that total CAS is individual to each patient
and organ offer, whereas the SRTR 2023 lung report uses WLAUC to characterize
medical urgency. Total CAS should not be substituted for a severity score:

- https://www.hrsa.gov/optn/patients/resources/lung/lung-allocation-faqs
- https://srtr.hrsa.gov/adr/2023/Lung/

## Implementation requirements

Obtain effective-dated calculated LAS and match/exception LAS separately, plus
post-transition calculated medical-urgency points or WLAUC and exception-adjusted
medical-urgency values. Request policy/formula versions and linkage IDs for every
registration, not only candidates who received offers. Include pre-LAS and
pediatric status applicability in the data specification.

Split each candidate's risk intervals at score-effective dates, pollution
boundaries and the March 9, 2023 transition. Use only information available at
each interval start; do not backfill from a future offer, transplant or score
update. Define same-day ordering explicitly when timestamps are absent. Carry
values forward only within their applicable score system, with missingness and
staleness audited. Use separate era-specific score coefficients/models and retain
the primary cohort's registration linkage rules. Do not pool LAS and total CAS
or present the existing baseline proxy as a time-updated official score.

Request candidate-level official LAS history before March 9, 2023, and CAS
medical-urgency component/history thereafter, with effective dates and linkage
identifiers. Confirm pediatric applicability and distinguish total CAS from its
medical-urgency component. Fit separate era-specific score terms or models;
do not pool LAS and CAS as a single interchangeable score.

The requested source audit is complete. Official-score replacement remains
blocked by the missing validated score-history source. No scores have been
invented or relabeled, and primary results have not been overwritten. The
baseline primary model has no organ score, so this limitation concerns the
score-adjusted analyses rather than the primary baseline model.
