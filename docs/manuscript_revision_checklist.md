# Manuscript revision checklist

Status: analyses completed, September 23, 2026. Items 1-6 and 8 are complete.
Item 7 is audited but official-score replacement requires a validated lung
score-history source. PH violations and sparse-covariate convergence warnings
are reported rather than assumed absent. Primary results are preserved;
new analyses are in `output/revision_2026_checklist`.

## Decisions confirmed with the investigator

- Competing-risk sensitivity: recent listing cohort, follow-up capped at 3 years.
- Relisting sensitivity: combine overlapping same-organ registrations within a
  continuous waitlist spell; separate spells after a gap.
- Baseline exposure: previous 365 days for PM2.5/O3, previous 12 complete months
  for NO2, retaining annual approximations before monthly coverage.

## 1. Simultaneous listings

- [x] Derive overlapping registration counts from registration start/end dates.
- [x] Distinguish concurrent listing count from lifetime registration count and
  distinct-center count. The existing episode_listing_count is not concurrency.
- [x] Fit baseline-exposure Cox sensitivity with concurrency updated as listings
  open and close; preserve the primary baseline exposure and covariates.
- [x] Also report the count at entry and its distribution. At first-ever listing
  this may be almost uniformly one and provide little adjustment.
- [x] Never use maximum future concurrency as a baseline covariate.
- [x] Report estimates separately from the time-varying pollution models.

## 2. Short-term competing-risk regression

- [x] Proposed recent cohort: listings during 2020-2022, conditional on verifying
  an administrative outcome endpoint through December 2025. This provides three
  years of potential observation without requiring complete observed follow-up.
- [x] Fit Fine-Gray regression for death/deterioration, with transplant/improvement
  competing; other exits remain censored under the existing outcome definition.
- [x] Cap follow-up at 3 years; use prelisting exposure and primary baseline
  demographic/SVI covariates. Report subdistribution HRs, not cause-specific HRs.
- [x] Explicitly document center handling. Censoring groups alone are not a
  substitute for center-stratified subdistribution baseline hazards.
- [x] Use the entire eligible recent cohort, sequential organ fits, and runtime/
  memory checkpoints. Do not select a random subset to make the fit feasible.

## 3. Separate sequential spells

- [x] Confirmed current code assigns person_organ_episode = 1 to every same-organ
  registration for a person (code/50, registration_cohort construction).
- [x] Replace that assignment only in the sensitivity cohort with overlap-based
  connected spells, using a running maximum end date to handle chained overlaps.
- [x] Keep same-day touching registrations together unless date auditing supports
  a different boundary rule; report that convention.
- [x] Reset time zero, center, address, baseline covariates, and prelisting exposure
  for each spell. Preserve the current within-spell outcome-resolution rule.
- [x] Refit the primary adjusted baseline Cox models; report people, spells,
  relistings, events, and changes from the reference cohort.
- [x] Evaluate person-clustered standard errors for repeated spells as a distinct
  variance sensitivity; do not silently change the reference variance estimator.

## 4. AJ numbers at risk

- [x] Existing AJ outputs contain risk-set counts, but primary plots do not show
  aligned risk tables (code/78 and code/figures_tables/83).
- [x] Compute exact risk sets at each displayed axis tick, including time zero,
  for each organ, pollutant, quartile, and horizon (1, 3, 5, 10 years).
- [x] Put each table directly below its curve panel. Include all quartiles and
  increase figure height; test labels at both endpoints.
- [x] Define at risk as still observed and free of either event immediately before
  the time point. Do not use the Fine-Gray augmented risk set for AJ tables.
- [x] Do not extend a last observed risk-set count past the end of follow-up.

## 5. Policy eras

- [x] Prespecify organ-specific listing-era cut points and pollutant-by-era
  interaction tests; report Ns/events and CIs, not just within-era significance.
- [x] Heart: October 18, 2018 allocation change; distinguish adult policy relevance
  from pediatric candidates.
- [x] Kidney: December 4, 2014 KAS and March 15, 2021 distribution change.
- [x] Liver: February 4, 2020 acuity circles and July 13, 2023 urgency-score changes;
  review whether earlier Share 35/MELD-Na transitions merit additional eras.
- [x] Lung: March 9, 2023 LAS-to-CAS transition; review earlier LAS implementation
  and geographic changes before finalizing the complete period definitions.
- [x] Label these as listing-era analyses. Analyses of policy operating during
  follow-up instead need intervals split at the policy date for existing candidates.
- [x] Avoid too many underpowered eras, particularly the short post-CAS period.

## 6. Proportional-hazards diagnostics

- [x] Fit/save scaled Schoenfeld residual diagnostics for pollutant terms and
  other nonstratified covariates, with global tests for each organ/model.
- [x] Evaluate baseline and time-varying models separately. A time-varying exposure
  does not itself allow its coefficient to vary with time since listing.
- [x] Plot residual trends and coefficient-versus-time curves with uncertainty.
- [x] If meaningful nonproportionality is present, report prespecified time-band
  effects or exposure-by-follow-up-time interactions instead of claiming PH holds.
- [x] Interpret tests alongside effect magnitude and plots in these large cohorts.

## 7. LAS/CAS audit and correction

- [x] NOT currently separated in the model code: code/50 constructs a weighted
  lung_las_cas_component_score, a proxy rather than official LAS or CAS points.
- [x] code/79 updates heart/liver scores, but lung retains baseline organ_score.
  It therefore does not currently implement a longitudinal official lung score.
- [x] Table 1 code labels the lung measure as a proxy linear predictor.
- [ ] Locate and validate official LAS, CAS, and medical-urgency component sources
  and timestamps in SAF/supplemental data. Absence from the current extract is not
  proof that a field is unavailable in the source data.
- [x] Do not pool official LAS and CAS under one continuous coefficient or relabel
  the proxy.
- [ ] Specify era-specific official score terms/models and validate pediatric use
  after receiving a validated score-history source.
- [ ] Update score-adjusted models and descriptions only after this source audit;
  the primary baseline Model 2 has no organ score and is unaffected by that change.

## 8. ACS proxy versus established SVI

- [x] Obtain CDC/ATSDR 2022 national ZCTA-level SVI and its documentation.
- [x] Match to the proxy calculated from the corresponding observed ACS vintage,
  not a carried-forward/backward value; verify geography, components, and direction.
- [x] Report matched/unmatched ZCTAs, missingness, Spearman correlation, and a
  scatter/hexbin plot against overall SVI and its socioeconomic theme.
- [x] Report quartile agreement, weighted kappa, and highest-vulnerability-quartile
  concordance. Use one row per ZCTA; optionally add population/candidate-weighted
  results, explicitly identified as secondary analyses.
- [x] Treat this as convergent validity, not proof of interchangeability. Report
  weak agreement if found; do not tune the proxy to force a positive validation.

## Recommended order

1. Resolve lung score provenance and reconstruct registration timelines.
2. Add AJ risk tables and run PH diagnostics.
3. Run concurrency, separate-spell, and recent-cohort Fine-Gray sensitivities.
4. Run prespecified policy-era analyses and SVI comparison.
5. Assemble one supplemental sensitivity table and update manuscript wording.

## Primary external sources

- Heart allocation: https://www.hrsa.gov/optn/news-events/news/four-year-monitoring-report-adult-heart-allocation-now-available
- Kidney KAS: https://optn.transplant.hrsa.gov/news/revised-national-kidney-transplant-allocation-system-is-now-in-place
- Kidney distribution: https://optn.transplant.hrsa.gov/news/updated-monitoring-report-available-for-kidney-allocation-policies/
- Liver acuity circles: https://www.hrsa.gov/optn/policies-bylaws/policy-issues/liver-intestine-policy
- Liver urgency scores: https://www.hrsa.gov/optn/news-events/news/updates-medical-urgency-scoring-liver-transplant-candidates-effect
- Lung continuous distribution: https://www.hrsa.gov/optn/professionals/resources/heart-lung/lung-continuous-distribution-policy
- CDC SVI ZCTA availability: https://www.atsdr.cdc.gov/place-health/php/svi/svi-news-updates.html
- CDC SVI documentation: https://www.atsdr.cdc.gov/place-health/php/svi/svi-data-documentation-download.html

## Execution findings

- Twelve recent-cohort Fine-Gray fits completed without subsampling, alongside
  twelve matched recent-cohort cause-specific Cox comparators. Both are attenuated
  relative to the full-period thoracic models; do not describe all sensitivities
  as confirming the full-period estimates.
- Separate-spell and concurrent-listing models preserve higher PM2.5/NO2 hazards
  in heart and lung. Model-based and person-clustered spell variances are separate.
- AJ figures include 1-, 3-, 5- and 10-year risk tables, with PM2.5/NO2 together and
  O3 separately. Endpoint text clipping was corrected during visual review.
- Policy eras additionally include MELD-Na (2016-01-11) and lung geographic
  allocation (2017-11-24). Interaction tests and counts accompany each estimate.
- The observed 2022 ACS source vintage was verified. Among 32,428 ZCTAs with valid
  official ranks, Spearman correlations were 0.767 with overall SVI and 0.858 with
  its socioeconomic theme. Component missingness and quartile agreement are saved.
- Lung audit details and the outstanding data requirement are in
  `docs/lung_score_source_audit.md`. PTR match-specific scores were not relabeled
  as longitudinal medical-urgency scores.
- All 24 baseline/time-varying primary estimates and sample/event counts were
  reproduced. All 48 AJ time-zero risk counts match the original curve cohorts.
- All 34 final PDFs were rendered for visual review; paired high-resolution PNGs
  passed nonblank checks. Diagnostic plots separate raw residuals from smoothed
  coefficient trends so uncertainty is legible.
- Supplemental CSV/TSV/HTML tables include convergence-warning flags and an audit
  table. Warnings involving sparse nuisance categories are not silently suppressed.
