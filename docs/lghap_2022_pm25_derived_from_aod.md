# LGHAP 2022 Model-Derived PM2.5 From AOD

LGHAP v2 currently provides native daily global PM2.5 grids through 2021. For
2022, the local archive contains LGHAP daily AOD grids, not native PM2.5 grids.
The 2022 PM2.5 workflow therefore produces a clearly labeled model-derived
extension rather than a native LGHAP PM2.5 release.

## Scripts

- `code/pollution_aggregation/45_aggregate_lghap_daily_aod_to_zcta.R` aggregates LGHAP daily AOD to
  2020 ZCTA5 polygons using the same rasterized ZCTA workflow and nearest-cell
  fill audit used for the native daily PM2.5 aggregation.
- `code/pollution_aggregation/46_model_derive_lghap_2022_pm25_from_aod.R` calibrates a ZCTA-level
  AOD-to-PM2.5 model using paired 2021 LGHAP AOD and native 2021 LGHAP PM2.5,
  plus 2005-2021 ZCTA monthly PM2.5 climatology, then predicts 2022 PM2.5 from
  2022 AOD.

## Required Inputs

- Native daily PM2.5 ZCTA aggregates for 2005-2021:
  `data/processed/lghap_pm25_zcta_daily/lghap_pm25_zcta_daily_YYYY_MM.csv.gz`
- LGHAP 2021 daily AOD monthly archives from Zenodo record `8301379`, or
  pre-aggregated 2021 AOD ZCTA files.
- Local 2022 LGHAP AOD archive:
  `C:/Users/Peter Graffy/Downloads/12702534.zip`

## Example Run

Aggregate 2022 AOD from the local outer archive:

```powershell
$env:LGHAP_DAILY_AOD_YEAR='2022'
$env:LGHAP_DAILY_AOD_OUTER_ZIP='C:/Users/Peter Graffy/Downloads/12702534.zip'
Rscript code\45_aggregate_lghap_daily_aod_to_zcta.R
```

Aggregate 2021 AOD for calibration from Zenodo:

```powershell
$env:LGHAP_DAILY_AOD_YEAR='2021'
$env:LGHAP_DAILY_AOD_ZENODO_RECORD='8301379'
Rscript code\45_aggregate_lghap_daily_aod_to_zcta.R
```

Fit the calibration model and predict 2022 PM2.5:

```powershell
Rscript code\46_model_derive_lghap_2022_pm25_from_aod.R
```

## Outputs

Daily model-derived 2022 PM2.5 ZCTA rows are written one month per file:

`data/processed/lghap_pm25_zcta_daily_derived/lghap_pm25_zcta_daily_derived_from_aod_2022_MM.csv.gz`

QA and model files:

- `data/processed/lghap_pm25_zcta_daily_derived/lghap_pm25_from_aod_calibration_model_2021.rds`
- `data/processed/lghap_pm25_zcta_daily_derived/lghap_pm25_from_aod_calibration_coefficients_2021.csv`
- `data/processed/lghap_pm25_zcta_daily_derived/lghap_pm25_from_aod_calibration_qa_2021.csv`
- `data/processed/lghap_pm25_zcta_daily_derived/lghap_pm25_zcta_daily_derived_from_aod_qc_2022_all_months.csv`
- `output/figures/lghap_pm25_derived_2022/lghap_pm25_from_aod_calibration_scatter_2021.png`

## Citation

Bai, K., Li, K., Shao, L., Li, X., Liu, C., Li, Z., Ma, M., Han, D., Sun, Y.,
Zheng, Z., Li, R., Chang, N.-B., and Guo, J.: LGHAP v2: a global gap-free
aerosol optical depth and PM2.5 concentration dataset since 2000 derived via big
Earth data analytics, Earth Syst. Sci. Data, 16, 2425-2448,
https://doi.org/10.5194/essd-16-2425-2024, 2024.
