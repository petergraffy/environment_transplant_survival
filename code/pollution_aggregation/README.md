# Pollution aggregation scripts

This folder contains scripts used to build, QA, package, and upload ZCTA-level
air pollution exposure surfaces. These scripts are upstream of the transplant
waitlist survival analyses in the parent `code/` directory.

Main workflows:

- Historical annual/monthly PM2.5, O3, and NO2 ZCTA aggregation:
  `aggregate_pollution_to_zcta.R`, `fill_missing_zip_pollution_values.R`,
  `25_package_air_pollution_release_assets.R`.
- LGHAP PM2.5 and AOD-derived PM2.5 ZCTA products:
  `39_aggregate_lghap_pm25_to_zcta.R`,
  `40_aggregate_lghap_daily_pm25_to_zcta.R`,
  `45_aggregate_lghap_daily_aod_to_zcta.R`,
  `46_model_derive_lghap_2022_pm25_from_aod.R`.
- Daily ozone ZCTA products:
  `52_aggregate_daily_o3_to_zcta.R`,
  `54_package_daily_o3_release_assets.R`.
- Monthly NO2 ZCTA products:
  `56_aggregate_monthly_no2_to_zcta.R`,
  `57_package_monthly_no2_release_assets.R`.
- QA maps and raw-grid inspections:
  `inspect_pollution_nc.R`, `plot_raw_pollution_grids.R`,
  `make_pollution_zip_maps.R`, and related map scripts.

Run scripts from the repository root so relative paths such as `data/`,
`output/`, and `code/r_runtime.R` resolve correctly.
