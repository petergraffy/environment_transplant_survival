# Air Pollution ZCTA Release Assets

Final filled air-pollution exposure surfaces are archived as GitHub Release
assets:

`https://github.com/petergraffy/environment_transplant_survival/releases/tag/air-pollution-zcta-v1`

## Assets

- `air_pollution_zcta_pm25_monthly_2005_2023.parquet`: monthly PM2.5 by 2020
  ZCTA, value column `pm25_ug_m3`.
- `air_pollution_zcta_o3_monthly_2005_2023.parquet`: monthly ozone by 2020
  ZCTA, value column `o3_ppb`.
- `air_pollution_zcta_no2_annual_2005_2024.parquet`: annual NO2 by 2020
  ZCTA, value column `no2`.
- `air_pollution_zcta_parquet_manifest.csv`: row counts, coverage,
  completeness, sizes, and SHA-256 checksums.
- `fill_summary.csv` and `fill_audit_summary.csv`: nearest-fill audit metadata.

All three Parquet files include `zip`, `pollutant`, `year`, `month`,
`temporal_resolution`, the pollutant value column, `value_source`, and
`fill_distance_km`.

PM2.5 and ozone currently cover 2005-2023 in the local source files; NO2 covers
2005-2024.
