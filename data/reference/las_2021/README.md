# Final pre-CAS LAS reference

The 365-day baseline survival tables are the **new** tables (PDF pages 17-18 and 21-22, one-based) in the official OPTN 2021 implementation notice:
https://www.hrsa.gov/sites/default/files/hrsa/optn/updated-cohort-for-calculation-of-the-las.pdf

Coefficients and diagnosis refinements are verified against the July 2021 Executive Committee notice, implemented September 30, 2021:
https://www.hrsa.gov/sites/default/files/hrsa/optn/20210730_executive_committee_summary.pdf

The research sensitivity applies this fixed formula retrospectively, not historical allocation-policy versions. It is not a recovered official LAS. Missing inputs use the published component-specific defaults. No future transplant labs or undated most-recent creatinine are used. Unmapped diagnoses and children younger than 12 have unavailable LAS, represented by a separate model indicator. Unmapped codes are reported, not silently assigned a diagnosis group. Sarcoidosis with missing mean pulmonary pressure is group A per policy.

Carry-forward beyond the policy expiration window is a deliberate user-requested research approximation. No unobserved clinical change is invented. CAS is dated offer-level PTR_TOT_SCORE, not a candidate-only medical urgency score.
