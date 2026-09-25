# LGHAP Daily PM2.5 ZCTA Parquet Release Assets

These files contain daily LGHAP PM2.5 exposures aggregated to 2020 ZCTA5 polygons for CONUS.
Coverage: 2005-2024.

Assets:

- `lghap_pm25_zcta_daily_YYYY.parquet`: yearly daily ZCTA PM2.5 values, one row per ZCTA-date.
- `lghap_pm25_zcta_daily_parquet_manifest.csv`: row counts, date coverage, completeness, file sizes, and SHA-256 checksums.
- `lghap_pm25_zcta_daily_qc_all_years.csv`: monthly completeness and distribution QC from the source aggregation.

Columns:

- `zip`
- `date`
- `pollutant`
- `year`
- `month`
- `temporal_resolution`
- `pm25_ug_m3`
- `value_source`
- `fill_distance_m`
- `fill_cell`

The `pm25_ug_m3` column is daily PM2.5 in micrograms per cubic meter. `value_source`, `fill_distance_m`, and `fill_cell` provide audit metadata for nearest-cell fills. All packaged years have zero missing ZCTA-day PM2.5 values. Join using `zip` (a five-character string) and `date`.

The 2022-2024 additions were aggregated from downloaded LGHAP daily PM2.5 NetCDF files using the same pipeline and 2020 ZCTA boundaries. The 2022 asset contains LGHAP PM2.5, not the separate locally derived AOD-to-PM2.5 estimates. Leap days are included.

Citation:

Bai, K., Li, K., Shao, L., Li, X., Liu, C., Li, Z., Ma, M., Han, D., Sun, Y., Zheng, Z., Li, R., Chang, N.-B., and Guo, J.: LGHAP v2: a global gap-free aerosol optical depth and PM2.5 concentration dataset since 2000 derived via big Earth data analytics, Earth Syst. Sci. Data, 16, 2425-2448, https://doi.org/10.5194/essd-16-2425-2024, 2024.

Also cite this repository/release for the ZCTA aggregation workflow. The reproducible aggregation script is `code/pollution_aggregation/40_aggregate_lghap_daily_pm25_to_zcta.R` and this packaging script is `code/pollution_aggregation/41_package_lghap_daily_pm25_release_assets.R`.
