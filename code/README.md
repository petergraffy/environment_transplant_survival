# Current primary analysis: lung transplant RMST benefit

Run all scripts from the repository root. The primary endpoint is the difference
in expected days alive over 365 days with versus without lung transplantation,
across PM2.5, NO2 and O3 exposure levels.

```powershell
Rscript code/run_primary_lung_analysis.R
```

Use `--reuse-fitted` to regenerate RMST curves and confidence intervals from
existing fitted models and matching bootstrap checkpoints. Without that flag,
the runner verifies/prepares cohort inputs and refits the six models.

| Script | Role |
|---|---|
| `00_prepare_lung_primary_inputs.R` | Prepare/verify established private cohort inputs |
| `91_lung_pollution_survival_benefit.R` | Pre/posttransplant Cox models and supporting survival probabilities |
| `93_lung_pollution_rmst.R` | Primary one-year RMST point estimates and figures |
| `94_lung_pollution_rmst_bootstrap.R` | Joint person bootstrap, RMST confidence intervals and figures |
| `92_lung_pollution_benefit_bootstrap.R` | Optional secondary survival-probability bootstrap |
| `lung_survival_benefit_helpers.R` | Shared model, prediction and exact RMST integration functions |

`r_runtime.R`, `saf_paths.R`, the daily/rolling pollution readers, exposure
preparation scripts and `06_build_community_covariates.R` remain shared inputs.
`pollution_aggregation/` contains upstream release construction.

The older analysis code is in [`archive/pre_lung_rmst/`](archive/pre_lung_rmst/README.md).
Its relative layout is preserved, and executable references have been updated.
The current cohort adapter reuses its established cohort-building functions;
the archived adverse-outcome models are not the primary analysis.

See [full methods and validation instructions](../docs/lung_pollution_survival_benefit.md)
and [primary results](../output/lung_pollution_survival_benefit/README.md).
