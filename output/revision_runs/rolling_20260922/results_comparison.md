# Rolling exposure revision results

Baseline PM2.5 and O3 now use the 365 days before listing. NO2 uses the 12 complete calendar months before the listing month, day weighted, with annual approximations before 2019. Daily PM2.5 coverage also extends from 2021 to 2024. These comparisons therefore combine changes in exposure definition and eligibility; they do not isolate one change at a time.

## baseline primary Cox models

| Organ | Pollutant | Previous HR | Updated HR (95% CI) | Previous N | Updated N |
|---|---|---:|---:|---:|---:|
| HR | pm25 | 1.253 | 1.259 (1.196-1.326) | 63,342 | 74,028 |
| HR | o3 | 1.007 | 0.992 (0.926-1.064) | 79,381 | 74,028 |
| HR | no2 | 1.151 | 1.175 (1.130-1.223) | 79,381 | 79,381 |
| KI | pm25 | 0.922 | 0.929 (0.916-0.942) | 492,197 | 564,179 |
| KI | o3 | 1.051 | 1.051 (1.033-1.068) | 600,141 | 564,179 |
| KI | no2 | 0.912 | 0.904 (0.895-0.915) | 600,141 | 600,141 |
| LI | pm25 | 1.018 | 1.025 (1.002-1.049) | 181,559 | 208,592 |
| LI | o3 | 1.023 | 1.034 (1.005-1.064) | 222,559 | 208,592 |
| LI | no2 | 0.986 | 0.997 (0.979-1.015) | 222,559 | 222,559 |
| LU | pm25 | 1.098 | 1.150 (1.076-1.229) | 40,175 | 46,788 |
| LU | o3 | 0.958 | 0.970 (0.892-1.056) | 50,401 | 46,788 |
| LU | no2 | 1.117 | 1.132 (1.076-1.190) | 50,401 | 50,401 |

## timevarying primary Cox models

| Organ | Pollutant | Previous HR | Updated HR (95% CI) | Previous N | Updated N |
|---|---|---:|---:|---:|---:|
| HR | pm25 | 1.140 | 1.144 (1.107-1.181) | 61,449 | 76,760 |
| KI | pm25 | 0.994 | 0.985 (0.976-0.995) | 486,543 | 592,129 |
| LI | pm25 | 1.032 | 1.037 (1.021-1.053) | 179,632 | 218,837 |
| LU | pm25 | 1.132 | 1.142 (1.095-1.191) | 38,786 | 48,311 |
| HR | o3 | 1.004 | 1.004 (0.983-1.026) | 76,760 | 76,760 |
| KI | o3 | 0.986 | 0.986 (0.980-0.992) | 592,129 | 592,129 |
| LI | o3 | 0.988 | 0.988 (0.978-0.998) | 218,837 | 218,837 |
| LU | o3 | 0.965 | 0.965 (0.938-0.992) | 48,311 | 48,311 |
| HR | no2 | 1.284 | 1.284 (1.234-1.336) | 82,122 | 82,122 |
| KI | no2 | 0.897 | 0.897 (0.886-0.908) | 628,132 | 628,132 |
| LI | no2 | 1.014 | 1.014 (0.996-1.032) | 232,844 | 232,844 |
| LU | no2 | 1.350 | 1.350 (1.283-1.422) | 51,932 | 51,932 |

## AJ cumulative incidence of death/deterioration

| Organ | Years | Pollutant | Quartile | Previous (%) | Updated (%) |
|---|---:|---|---|---:|---:|
| Heart | 1 | NO2 | Q1 lowest | 8.50 | 7.89 |
| Heart | 1 | NO2 | Q4 highest | 10.59 | 10.97 |
| Heart | 1 | PM2.5 | Q1 lowest | 8.62 | 8.06 |
| Heart | 1 | PM2.5 | Q4 highest | 11.83 | 11.65 |
| Kidney | 10 | NO2 | Q1 lowest | 26.22 | 25.37 |
| Kidney | 10 | NO2 | Q4 highest | 28.43 | 28.88 |
| Kidney | 10 | PM2.5 | Q1 lowest | 24.47 | 24.09 |
| Kidney | 10 | PM2.5 | Q4 highest | 29.08 | 28.78 |
| Liver | 3 | NO2 | Q1 lowest | 16.02 | 15.15 |
| Liver | 3 | NO2 | Q4 highest | 21.33 | 21.78 |
| Liver | 3 | PM2.5 | Q1 lowest | 18.62 | 17.40 |
| Liver | 3 | PM2.5 | Q4 highest | 20.68 | 20.65 |
| Lung | 3 | NO2 | Q1 lowest | 8.56 | 8.13 |
| Lung | 3 | NO2 | Q4 highest | 15.24 | 15.76 |
| Lung | 3 | PM2.5 | Q1 lowest | 11.08 | 9.80 |
| Lung | 3 | PM2.5 | Q4 highest | 14.79 | 14.14 |

## Baseline two-pollutant thoracic models

| Organ | Term | Previous HR | Updated HR (95% CI) | Updated P |
|---|---|---:|---:|---:|
| HR | pm25_prior_5ug | 1.232 | 1.191 (1.121-1.267) | 2.012e-08 |
| HR | no2_prior_10ppb | 1.027 | 1.084 (1.034-1.136) | 0.0008779 |
| LU | pm25_prior_5ug | 1.066 | 1.092 (1.010-1.181) | 0.02648 |
| LU | no2_prior_10ppb | 1.045 | 1.078 (1.016-1.144) | 0.01357 |

## Scope and interpretation

PM2.5 HRs are per 5 ug/m3; NO2 and O3 HRs are per 10 ppb.

Scope: current baseline and monthly time-varying models, all three baseline adjustment tiers, subgroup and sensitivity analyses, severity associations, AJ curves, and their current figures/tables were regenerated. Superseded whole-waitlist-average, KM, and Fine-Gray analyses remain archived rather than being relabeled as current results.

Table 1: candidate counts, demographics, clinical characteristics, and outcome rows are unchanged. The three exposure rows now describe prelisting exposures rather than whole-waitlist averages.

The very large relative changes in the full-coefficient audit concern the kidney missing-race indicator, whose previous HR was near zero; these are nuisance-category coefficients, not pollutant associations. All pollutant estimates are separately compared in the baseline, model-tier, time-varying, subgroup, and sensitivity files.

Figure QA: 48 updated PNG exports passed image checks and 48 updated PDFs opened successfully. Contact sheets are in figure_qa/. The main forest PDF was rendered separately for visual inspection. The full AJ and LVAD/dialysis AJ sets are in rolling_365d_20260922 subfolders to avoid locked previous exports.

Interpretation caution: the existing map calculation uses 1 - exp(-H_death(365)) from a cause-specific Cox model. This is a net-risk transformation, not the cumulative incidence of death/deterioration before transplant. The latter requires the competing-event hazard as well. This exposure revision preserves the map calculation; it does not resolve that separate estimand issue.

All model tiers, full coefficients, subgroups, sensitivities, severity associations, and AJ horizons have separate comparison CSVs in this directory. A changed significance classification alone does not establish a statistically significant difference between the old and new estimates.
