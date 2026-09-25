# Exploratory lung LAS/CAS sensitivity

## Scope

This is an additional sensitivity, not a replacement for the primary analysis. The original models, figures and Table 1 are unchanged. All single-pollutant, two-/three-pollutant and age/sex/race subgroup time-varying lung models were fitted. No-center and multiorgan sensitivity variants were also fitted: 41 models total, providing 44 pollutant coefficients.

## Sources and reconstruction

The fixed September 2021 LAS equation (the final pre-CAS revision) was applied retrospectively. The reference curves and official sources are in `data/reference/las_2021/`. This is not an era-specific reconstruction of the official allocation score. Its equation is:

`LAS = 100 * (PTAUC - 2 * WLAUC + 730) / 1095`

Each AUC is the sum of days 0 through 364 of the published baseline survival probability raised to `exp(linear predictor)`. Coefficients, diagnosis mapping and transformations are explicit in `code/revision_2026/las_2021.R`. This is not sample min-max scaling.

Candidate fields: CAN_AGE_AT_LISTING, CAN_DGN, CAN_BMI, CAN_HGT_CM, CAN_WGT_KG, CAN_FUNCTN_STAT, CAN_VENTILATOR (fallback CAN_ON_VENTILATOR), CAN_AT_REST_O2, CAN_PCO2, CAN_SIX_MIN_WALK, CAN_PULM_ART_SYST, CAN_PULM_ART_MEAN, CAN_CARDIAC_OUTPUT. Cardiac index is cardiac output divided by Mosteller BSA, `sqrt(height_cm * weight_kg / 3600)`. The ventilation flag approximates the policy's hospitalized continuous ventilation term; hospitalization and ventilation timing cannot be fully verified. Functional independence is approximated by old code 1 or Karnofsky/Lansky >=70. New/unmapped diagnosis codes (including 1616, 1617 and 1620, absent from the supplied format catalog) are reported rather than guessed. Candidates younger than 12 or with an unmapped diagnosis have unavailable LAS, modeled with an indicator, rather than being dropped.

No temporally verifiable LU creatinine or bilirubin measurement was available. CAN_MOST_RECENT_CREAT has no test timestamp and was not assigned backward to listing. Transplant-recipient laboratory measurements were not used. The official defaults were retained: creatinine 0.1 mg/dL in the waiting-list component and 40 mg/dL in the adult post-transplant component (neither contributes under age 18); bilirubin 0.7 mg/dL. Other component-specific defaults are shown in the Table 1 addendum. Serial PCO2 increase was not inferred without dated tests. FVC, FEV1, diabetes, ECMO and corticosteroids are not additional terms in the final LAS formula used here.

Inputs are assumed to describe status on their candidate registration date, as in the existing candidate-field analysis. At a later registration, observed inputs replace prior inputs and missing inputs carry forward; if no clinical input is supplied, the score remains unchanged. Same-day candidate records prioritize the baseline registration, then the most complete form, then PX_ID. Updates after initial listing become effective the next day. Thus no later observation is applied before its availability. Carry-forward is indefinite at the user's request, unlike the official policy's expiration rules. Baseline score is retained if there is no subsequent usable record. No unobserved laboratory trajectories are interpolated.

Beginning March 9, 2023, the post-transition score is PTR_TOT_SCORE from lung offer records, linked across candidate registrations using PX_ID and PERS_ID. Overlapping annual files are pooled. Each candidate-day uses the median of distinct observed score values; it becomes effective the following day and is carried forward. No backward filling from future offers is allowed. This is offer-level total CAS, not isolated candidate medical urgency, and includes donor-dependent allocation information. Before the first observed CAS, the post-transition score is unavailable and a missingness indicator is used.

## Models

For each pollutant, monthly exposure is unchanged from the original time-varying analysis: daily-derived monthly PM2.5/O3 through December 2024 and monthly NO2 with annual fallback through December 2025. Monthly intervals are additionally split at score updates and the policy transition. A candidate-day appears once; events remain on the terminal interval. The primary sensitivity formula is:

`Surv(tstart, tstop, adverse_event) ~ exposure_terms + age + sex + race + ACS_vulnerability_proxy + post_CAS + LAS_pretransition/10 + CAS_posttransition/10 + LAS_unavailable + CAS_unavailable + strata(first_listing_center)`

Inactive score terms are zero outside their respective eras; an unavailable active score is zero with its separate unavailable indicator. The old organ score is not included alongside these terms. Listing year is not included. Efron ties and model-based standard errors are used. A constant nuisance covariate is omitted within a subgroup. Age-subgroup models omit age; sex/race subgroup models use attained age and omit the stratifying covariate, matching existing subgroup conventions. Minimum 50 events and 2 centers are required. The three under-18 subgroup fits produce an infinite-coefficient warning for a nuisance term and are flagged, not treated as stable inferential results.

## Validation and interpretability

Unit tests cover diagnosis categories, coefficient arithmetic, missing-value substitutions, the AUC normalization, and under-12 applicability. Interval checks confirm contiguous, nonoverlapping time and preservation of person-time/events. Original-score models refitted on the split dataset reproduce the existing HRs to within 0.000001 with identical candidate/event counts. Separate original-score models with the transition indicator quantify the effect of era adjustment alone.

There are 52,613 lung candidates; 50,873 have a computable reconstructed baseline LAS, and 2,961 have a usable subsequent pre-transition listing score record. The reconstructed baseline LAS median is 0.39 (IQR, 0.24-1.95). Oxygen, walk distance and PCO2 are missing for approximately two-thirds or more of the cohort, and dated creatinine/bilirubin are wholly unavailable. Consequently this estimate is heavily driven by least-beneficial policy defaults and should not be described as clinically validated LAS. Offer-driven CAS availability can also be informative. Attenuation of pollution estimates in this sensitivity is not evidence of mediation or proof that full adjustment for clinical severity removes the association.

## Reproduction and outputs

Run with R from the repository root:

1. `code/revision_2026/19_inspect_las_inputs.R`
2. `code/tests/test_las_2021.R`
3. `code/revision_2026/20_lung_las_cas_sensitivity.R`
4. `code/revision_2026/21_lung_las_cas_report.R`

Private derived data remain in the ignored `output/revision_2026_checklist/cache/` directory. Aggregate reports are in `output/revision_2026_checklist/lung_las_cas_sensitivity/`: publication-oriented HTML, pollutant results, full log-coefficients and standard errors, model formulas/warnings, input Table 1 addendum, categorical distributions, score coverage/update counts, diagnosis audit and validation benchmarks. Interval/CAS caches are reusable only for this version of the reconstruction; delete those specific derived caches before rerunning after changing score logic. Do not delete source SAF data.
