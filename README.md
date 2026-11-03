# Air pollution and lung transplant survival benefit

The primary analysis estimates how **one-year lung transplant survival benefit,
expressed as restricted mean survival time (RMST) gained**, varies with PM2.5,
NO2 and O3 exposure. It uses separate pretransplant and posttransplant
time-varying Cox mortality models and paired predictions over the same recipients.

The earlier four-organ waitlist death/deterioration analyses are preserved in
[`code/archive/pre_lung_rmst/`](code/archive/pre_lung_rmst/README.md). Historical
outputs remain in their original locations for provenance.

## Primary workflow

Run from the repository root with the existing R environment, local Q1 2026
SRTR SAFs, pollution releases and ACS community covariates:

```powershell
Rscript code/run_primary_lung_analysis.R
```

1. `code/00_prepare_lung_primary_inputs.R`: prepare or verify the private,
   audited LAS-era cohort inputs, retaining the established cohort definitions.
2. `code/91_lung_pollution_survival_benefit.R`: fit six time-varying Cox models.
3. `code/93_lung_pollution_rmst.R`: integrate fitted survival over 365 days.
4. `code/94_lung_pollution_rmst_bootstrap.R`: add pointwise 95% intervals by
   jointly resampling people, refitting both phases, and reintegrating RMST.
5. Validate the model outputs and RMST confidence intervals.

To regenerate RMST figures and restore confidence intervals from the existing
fitted models and matching bootstrap checkpoints:

```powershell
Rscript code/run_primary_lung_analysis.R --reuse-fitted
```

Bootstrap defaults: 200 attempts, up to 12 local workers, seed 20261102.
`LUNG_RMST_BOOTSTRAPS` and `LUNG_RMST_WORKERS` override those settings.
Completed replicates are checkpointed, so interrupted runs can resume. A full
refit changes the model-cache signature and requires new bootstrap estimates.

## Study definition

- Adults listed for lung transplantation from January 1, 2006 through March 8,
  2023, using the established first-listing/deduplicated cohort and exclusions.
- All-cause mortality, including ascertainable deaths after delisting and
  perioperative deaths. Pretransplant follow-up ends at the LAS transition;
  postoperative follow-up can extend through December 31, 2024.
- Pollution updates from completed historical months; older NO2 uses a lagged
  annual fallback. Listing residence and community vulnerability remain fixed.
- Adjustment: age, sex, race, ZCTA vulnerability, official LAS with its missing
  indicator when applicable, and listing-center strata. LAS updates before
  transplantation and is fixed at transplantation afterward.
- RMST benefit is expected days alive after transplant minus conditional
  expected days alive without transplant during the same next 365 days.

These are adjusted observational model estimates. The models do not identify
causal treatment effects or effects of reducing pollution. Confidence intervals
are pointwise and conditional on the fixed recipient reference population.
One-year survival-probability differences are supporting outcomes; script 92
provides their separate, optional bootstrap workflow.

## Results

The current primary results are in
[`output/lung_pollution_survival_benefit/`](output/lung_pollution_survival_benefit/README.md):

- [RMST benefit with 95% intervals](output/lung_pollution_survival_benefit/rmst_1yr_benefit.png)
- [RMST with and without transplant](output/lung_pollution_survival_benefit/rmst_1yr_components.png)
- [Accumulated benefit through one year](output/lung_pollution_survival_benefit/rmst_1yr_time.png)
- [Exposure contrasts and confidence limits](output/lung_pollution_survival_benefit/rmst_1yr_exposure_contrasts.csv)

The November 2, 2026 run used 30,186 paired reference recipients. Of 200
bootstrap attempts, 174 succeeded for each pollutant; 25 lacked complete
follow-up support and one had an unstable coefficient. The published tables
and plots contain the corresponding 95% percentile intervals.

Full definitions, assumptions, individual commands, source requirements and
validation are in the [primary analysis methods](docs/lung_pollution_survival_benefit.md).

## Repository layout and restricted data

- `code/`: current lung analysis, shared runtime/path helpers and input preparation.
- `code/tests/`: current lung and shared exposure-input checks.
- `code/pollution_aggregation/`: upstream exposure construction and release tools.
- `code/archive/pre_lung_rmst/`: historical models, figures, revision scripts,
  exploratory work and their tests. Cohort preparation still reuses documented
  functions from this archive.
- `docs/`: current methods and source documentation; older analysis notes are
  historical rather than current primary specifications.
- `output/`: aggregate results, figures and historical output bundles.

SAFs and person-level derived data are not redistributed. Private cohorts,
model objects and bootstrap working inputs stay in ignored `cache/` folders.
Exposure-release Parquet files and other large local inputs are also ignored.
