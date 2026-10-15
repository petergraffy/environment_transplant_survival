# New official lung score supplement: source audit

Source: Q1 2026 SAF/Lung_supplement2603. Audited October 15, 2026.
Primary models, tables, and figures have NOT been replaced by this audit.

## LAS

las_saf2603.sas7bdat has 13,019,510 records for 62,000 PX_ID values.
CANHX_LAS_SCORE_CALC is labeled calculated LAS; CANHX_LAS_SCORE_MATCH is
the match score and may be a requested score. CANHX_LAS_SCORE_BEGIN_DT
is explicitly the date when the score should start to be used.
Calculated LAS is missing on all records dated in 2024-2026. The 2005-2007
record-weighted calculated-score medians are zero; the meaning of zeros,
same-day updates, and baseline coverage require checking before fitting.
This file provides a longitudinal pre-CAS source, not a post-transition CAS history.

## CAS limitation

cas_2603.sas7bdat has 9,813 rows for 9,812 registrations (9,763 people).
Fields are pers_id, px_id, CAS, CAS_dt, and DistanceInNM. Dates span
2023-03-09 to 2026-02-28. Donor distance ranges from 0 to 2,920.47 nautical miles.
There is essentially one record per registration, not a serial score history.

Of 9,774 rows linked by PERS_ID/PX_ID to lung candidate records, 9,771
belong to registrations with removal code 4 (transplantation), and 3 to
code 21 (death during transplantation). CAS dates need not be transplant
dates: median timing in code-4 registrations is 8 days after listing and
7 days before removal. Do not infer the score's definition from timing alone.

Linkage to the collapsed primary cohort found an in-follow-up CAS record
for 9,149 transplant/improvement outcomes but only 5 adverse outcomes.
This is a broad PERS_ID coverage audit, not an approved modeling join.
Denominators in cas_coverage_by_primary_outcome.csv include the entire
2005-2025 cohort, not only the CAS era. Different final outcomes can arise
from the existing registration-collapse rules.

## Required clarification before primary replacement

Confirm with SRTR the construction and meaning of CAS and CAS_dt,
whether this is a recipient-selected or donor-specific extract, and obtain
dated candidate-level CAS/medical-urgency history for all waitlist outcomes.
Do not backfill a later CAS into earlier follow-up or treat this missingness
as an ordinary baseline missing-data problem. Keep LAS and CAS as separate
covariates/scales; do not splice them into one continuous score.

Reproduction: code/revision_2026/29_audit_official_lung_scores.R.
Aggregate audits are saved here; row-level source data remain in the ignored cache.
