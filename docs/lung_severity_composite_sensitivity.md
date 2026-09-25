# Study-defined lung severity composite sensitivity

## Purpose

Use the same clinical composite before and after March 9, 2023, without mixing LAS and CAS or including donor-dependent offer scores. Existing primary outputs and prior exploratory LAS/CAS outputs are preserved. This composite is neither a validated LAS nor a CAS, and its units are not allocation points.

## Fixed definition chosen before model fitting

The score is a reduced linear predictor using the well-recorded fields overlapping the **2021 LAS waiting-list mortality** equation. It excludes age because age is separately adjusted in the regression. It excludes oxygen, walk distance and PCO2 because roughly two-thirds or more of baseline observations are missing; dated bilirubin and creatinine are unavailable. Cardiac index is not added: although reasonably available, it is in the final LAS post-transplant component, not its waiting-list mortality component. No coefficients are estimated using this study's outcomes or pollution associations.

The five contributions are summed after fixed imputation:

| Domain | Transformation and weight |
|---|---|
| Diagnosis | A: 0; B: 1.26319338239175; C: 1.78024171092307; D: 1.51440083414275. Add 0.40107198445555 for bronchiectasis in A; add 0.2088684500011 for other pulmonary fibrosis in D; add 1.39885489102977 for sarcoidosis in A or subtract 0.64590852776042 for sarcoidosis in D. |
| Functional dependence | 0.59790409246653 for some/total assistance; 0 for independence. This is a constant-shifted version of the published protective independence coefficient. |
| Ventilation | 1.57618530736936 for ventilator support; 0 otherwise. |
| Low BMI | 0.10744133677215 * max(20 - BMI, 0). |
| Pulmonary artery systolic pressure | A: 0.55767046368853 * max(PASP - 40, 0) / 10. B/C/D: 0.1230478043299 * max(PASP, 20) / 10. |

The code-to-diagnosis-group mapping is shared with the documented reconstruction. An unmapped diagnosis is not guessed. Sarcoidosis with missing mean PA pressure is treated as unknown group here, not automatically assigned group A. Functional independence means old code 1 or Karnofsky/Lansky 70-100; old codes 2/3 or Karnofsky/Lansky 10-60 indicate dependence. Unknown status is missing, not independent. Ventilation uses CAN_VENTILATOR with CAN_ON_VENTILATOR as fallback; it approximates continuous ventilation because detailed timing is unavailable.

Missing contributions receive the observed median of that contribution among candidates initially listed **before** March 9, 2023. Those medians are frozen across both eras. Five domain-specific missingness indicators enter models alongside the score; constant indicators are omitted. Thus missing components are neither silently set to healthy nor omitted from a variable-length denominator. This remains a pragmatic missing-data approach, not multiple imputation or proof that missingness is ignorable.

Source: final 2021 OPTN LAS waiting-list coefficients, https://www.hrsa.gov/sites/default/files/hrsa/optn/20210730_executive_committee_summary.pdf . Published weights originally supported allocation scoring in candidates at least 12 years old; use of this unvalidated study composite in younger candidates/subgroups is exploratory.

## Longitudinal handling

The first registration supplies baseline fields. A later registration can update observed components from that date forward; unobserved components carry their previous observation forward. Same-day forms prioritize the baseline record, then completeness, then PX_ID. Subsequent updates take effect the following day to avoid same-day lookahead. Carry-forward does not imply that actual physiology remained unchanged. Dates on these candidate records are registration dates rather than serial laboratory collection timestamps. The data do not provide dense within-registration clinical histories.

Monthly pollution exposures remain unchanged. Follow-up intervals split at a change in composite value or missingness and at the policy transition. No PTR_TOT_SCORE or other offer-derived information is used.

## Models and eras

Pooled time-varying cause-specific models use pollution term(s), age, sex, race, ACS vulnerability proxy, composite, missingness indicators and listing-center strata. Listing year and an era indicator are not added to this pooled counterpart, matching the primary model's temporal specification. Subgroups follow existing age/sex/race adjustment conventions. The original organ score is replaced, not included alongside the composite.

Two clearly distinguished era analyses are exported:

1. **Listing-era cohorts**, matching the earlier era approach: before March 9, 2023 vs on/after March 9, 2023. Pre-transition listings may have post-transition follow-up. The older three-era lung analysis with the 2017 geography change remains untouched.
2. **Calendar-era risk intervals**: follow-up is split on March 9, 2023; candidates may contribute to both periods, but never duplicate person-time. This evaluates associations during each allocation-policy period.

The cutoff corresponds to adoption of CAS on March 9, 2023 (https://srtr.hrsa.gov/adr/2023/Lung/). Each analysis provides separate era fits and a formal 1-df likelihood-ratio interaction test comparing a common pollution coefficient with era-specific coefficients, using center-by-era strata and otherwise common nuisance coefficients. Separate era models allow their nuisance coefficients to vary and need not equal the estimates from the interaction model exactly. Interaction is not a test of a causal policy effect. Post-transition precision is limited, particularly for PM2.5/O3 ending in 2024. Model-based SEs and Efron ties match the existing analyses.

## Reproduce

Run `code/tests/test_lung_composite.R`, then `code/revision_2026/22_lung_composite_sensitivity.R`, then `code/revision_2026/23_lung_composite_report.R`, then `code/tests/test_lung_composite_outputs.R`. This uses the existing raw lung cache and unchanged monthly exposure cache from the preceding audit. Person-level derived files stay in ignored `output/revision_2026_checklist/cache/lung_composite/`. Aggregate tables and the report are in `output/revision_2026_checklist/lung_composite_sensitivity/`.

## Run findings and numerical limitations

Pooled single-pollutant HRs were 1.09 (1.04-1.14) for PM2.5 per 5 ug/m3, 1.17 (1.11-1.23) for NO2 per 10 ppb, and 0.97 (0.94-1.00) for O3 per 10 ppb. Listing-era interaction P values were .214, .428 and .804, respectively; calendar-era values were .115, .421 and .233. All six interaction models converged without warnings. Primary exposure models replicated after interval splitting.

Of 53 separate model attempts, the three Asian subgroup fits and the PM2.5/O3 post-transition listing-only fits encountered numerical overflow. Seven other fits produced infinite-coefficient warnings, associated with sparse subgroup/post-transition covariates. These are explicitly flagged in the exported results; they are not treated as reliable independent estimates. The main era table therefore uses the estimable pooled interaction models with common nuisance coefficients. This change in presentation does not discard any attempted fit or change the model specification to restore statistical significance.

All 52,613 lung candidates were retained in preparation; model-specific exposure/covariate restrictions match the original cohorts. Only 2,973 have a subsequent record changing the score or missingness status. Baseline domain availability is approximately 92%-100%. The composite remains an unvalidated research adjustment, with limited measurement frequency and possible residual confounding.
