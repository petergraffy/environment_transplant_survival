# Manuscript revision results

Primary results were preserved. New outputs are in this directory; executable
code is in `code/revision_2026`. Effect increments are 5 ug/m3 for PM2.5 and
10 ppb for NO2/O3. Baseline sensitivities retain age, sex, race, ACS vulnerability
proxy and center strata, without listing year or organ score.

Some follow-up-band, policy-era and recent-cohort models have sparse categorical
nuisance-coefficient warnings. Affected estimates are flagged [a] in the
supplemental tables and the warning audit is provided. In the baseline/recent
audit these involve missing or rare race categories, not the pollutant term.
Race categories were not silently combined to eliminate warnings. Interpret
flagged models cautiously; they should not be described as convergence-clean.
All 18 time-varying-run warnings were reproduced and localized to follow-up-band
fits, not the full-period primary fits. The baseline/era audit localized 30
warnings to band or era fits; each recent-cohort regression type had three
kidney missing-race warnings. The supplement contains these 54 warning records.

The first time-varying diagnostic process exited with a post-fit cleanup error
after saving all 12 models. All 24 primary reproductions subsequently passed
independent checks, and warning-producing fits were independently rerun with
matching estimates, candidate counts and events. Patient-level cache files
remain git-ignored.

## Registration sensitivities

Separate-spell reconstruction yielded 87,231 heart, 690,923 kidney, 252,217 liver
and 55,206 lung spells before exposure and complete-case exclusions. There were
3,466 heart, 44,807 kidney, 13,215 liver and 2,486 lung person-organ records with
more than one spell. Positive gaps define separate spells; same-day touching and
overlapping registrations remain combined. Spells reset the baseline date,
center, address, covariates and preceding-year exposure.

Among positive-duration primary timelines, 1,651 heart, 80,425 kidney, 10,271 liver
and 1,026 lung records had concurrent registrations at some point. Entry counts
are not a substitute for these longitudinal counts. The concurrency model also
includes a no-open-listing indicator for gaps retained in the primary timeline.
These counts precede pollutant-specific eligibility restrictions.

The thoracic associations persisted:

| Organ/pollutant | Primary baseline HR | Separate spells HR | Concurrency-adjusted HR |
|---|---|---|---|
| Heart PM2.5 | 1.26 (1.20-1.33) | 1.30 (1.24-1.37) | 1.31 (1.25-1.38) |
| Heart NO2 | 1.18 (1.13-1.22) | 1.21 (1.16-1.25) | 1.22 (1.17-1.27) |
| Lung PM2.5 | 1.15 (1.08-1.23) | 1.14 (1.07-1.22) | 1.13 (1.06-1.21) |
| Lung NO2 | 1.13 (1.08-1.19) | 1.13 (1.07-1.19) | 1.13 (1.08-1.19) |

Parentheses show 95% confidence intervals. Separate-spell clustered standard
errors yielded very similar intervals; both variance estimators are reported.
Full results, including kidney, liver and ozone, are in `registration_models`.

## Recent-cohort competing risks

### Updated sensitivity: 2015-2019 listings

The current supplemental tables use all eligible candidates listed January 1,
2015 through December 31, 2019, with follow-up capped at 3*365.25 days. This
excludes 2020 listings but does not exclude pandemic-era follow-up. The original
2020-2022 results below are retained for provenance, not used in the current
supplemental sensitivity columns. Script: `16_pre_covid_competing_risks.R`.

There were 20,504 heart, 145,305 kidney, 56,090 liver and 13,532 lung candidates,
with 2,531, 16,667, 9,926 and 1,436 adverse events, respectively. Fine-Gray and
matched cause-specific models have identical eligibility and event counts.
Prior-year exposure and adjustment remain age, sex, race, ACS vulnerability
proxy and center strata; listing year and organ score are not included.

| Organ | PM2.5 subdistribution HR (95% CI) | NO2 | O3 |
|---|---|---|---|
| Heart | 1.20 (1.03-1.41) | 1.06 (0.94-1.19) | 1.04 (0.88-1.23) |
| Kidney | 0.91 (0.85-0.96) | 0.85 (0.81-0.89) | 1.09 (1.03-1.16) |
| Liver | 1.02 (0.95-1.09) | 0.94 (0.89-0.99) | 1.03 (0.96-1.11) |
| Lung | 1.09 (0.92-1.30) | 1.24 (1.07-1.42) | 0.96 (0.79-1.17) |

All three kidney models of each regression type had a possible infinite
missing-race coefficient (six warnings total). Pollutant coefficients are
finite; affected estimates remain flagged. Neither heart nor lung models had
convergence warnings. Changing the listing window is not a causal test of
COVID effects. Matched cause-specific heart PM2.5 was 1.14 (0.97-1.33) and
lung NO2 was 1.21 (1.05-1.40).

Full results and coefficient files: `finegray_2015_2019/`. Formatted table:
`supplement/competing_risks_2015_2019.html`. The original cohort was not altered.

### Original sensitivity: 2020-2022 listings (archived comparison)

All eligible 2020-2022 listings were included, with follow-up capped at three
years: 13,336 heart, 93,355 kidney, 36,056 liver and 8,238 lung candidates. This
is a deliberately different cohort from the full-period analyses, not a random
subsample. Transplant/improvement are competing events; other exits are censored.

Fine-Gray estimates were attenuated or inverse, rather than confirming the
full-period thoracic associations. Heart subdistribution HRs were 0.82
(0.64-1.05) for PM2.5 and 0.72 (0.52-0.99) for NO2; lung estimates were 0.91
(0.69-1.21) and 0.70 (0.47-1.03), respectively. Kidney NO2 was 0.86 (0.78-0.96).

Matched recent-cohort cause-specific estimates were also attenuated: heart
PM2.5 0.87 (0.67-1.14), heart NO2 0.73 (0.52-1.03), lung PM2.5 0.91
(0.67-1.24), and lung NO2 0.78 (0.53-1.15). Thus the difference from the
full-period findings is not attributable solely to the Fine-Gray estimand.
These estimates do not establish a protective causal effect of pollution.

## Listing eras and proportional hazards

Era interactions were nominally significant for heart PM2.5 (P=.039) and NO2
(P=.045), and all three kidney pollutants (PM2.5 P<.001, NO2 P<.001,
O3 P=.033). No interaction reached .05 for liver or lung. These exploratory
P values are not multiplicity-adjusted. Heart era analyses are adult-only;
post-CAS lung estimates have wide intervals.

All baseline global Schoenfeld tests reject proportional hazards. Pollutant
tests reject a constant coefficient for heart PM2.5/NO2, all kidney pollutants,
and liver PM2.5/NO2. Lung pollutant-specific tests do not reject at .05, although
nuisance-covariate violations remain. A nonsignificant test does not prove PH.

For illustration, heart baseline PM2.5 HRs were 1.35 (1.27-1.43) in year 0-1,
1.10 (0.97-1.25) in years 1-3, 1.05 (0.84-1.31) in years 3-5, and 0.91
(0.67-1.22) after year 5. The constant-coefficient HR should therefore be
described as a summary over follow-up, not a constant effect at every time.
Time-varying model diagnostics and analogous follow-up-band estimates are
reported separately in `timevarying_diagnostics`.

All 12 time-varying global tests also reject PH. Pollutant-specific tests reject
PH for heart PM2.5/NO2, kidney PM2.5/NO2, and all three liver pollutants, but not
for the lung pollutant terms. Heart time-varying NO2 HRs were 1.37 (1.32-1.44)
in year 0-1 and 1.10 (0.99-1.22) in years 1-3. These results should qualify the
constant-HR interpretation of the main models. All 12 baseline and all 12
time-varying primary point estimates and sample/event counts were reproduced.

## ACS vulnerability validation

The observed ACS source vintage was verified as 2022. Of 33,774 proxy ZCTAs,
33,642 matched official CDC geography; 32,428 had usable official SVI ranks.
Spearman correlation was 0.767 with overall SVI (bootstrap 95% CI,
0.763-0.773) and 0.858 with its socioeconomic theme (0.855-0.862).
Exact quartile agreement was 52.1% and 63.5%, respectively; linearly weighted
kappa was 0.556 and 0.681. This supports convergent validity, not equivalence.
Component missingness and the primary median-imputation approach are documented.
Bootstrap intervals do not account for spatial dependence.

## Figures and unresolved source requirement

`aj_with_risk_tables` contains 1-, 3-, 5- and 10-year figures with aligned
numbers at risk, PM2.5/NO2 together and ozone separately. Numbers represent
observed, event-free candidates immediately before each tick, not the augmented
Fine-Gray risk set.

The lung source audit found that official LAS and CAS are **not** currently
separated. Current score-adjusted models use a baseline component proxy; the
lung component is not longitudinally updated. A validated candidate-level LAS
and CAS medical-urgency history is required before replacing it. PTR match-total
scores were not substituted. See `docs/lung_score_source_audit.md`.

## Supplemental tables

The `supplement` directory contains CSV, TSV and HTML tables. HTML tables can be
opened in a browser and copied into Word. The wide sensitivity table puts the
primary, concurrency, separate-spell and recent-cohort estimates side by side;
the long version adds denominators and event counts. Additional tables cover
policy eras, PH tests, follow-up bands, registration counts and SVI validation.
