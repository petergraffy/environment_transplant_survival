# Daily Ozone ZCTA Parquet Release Assets

These files contain daily MDA8 ozone exposures aggregated to 2020 ZCTA5 polygons for CONUS.

Assets:

- `o3_zcta_daily_YYYY.parquet`: yearly daily ZCTA ozone values, one row per ZCTA-date.
- `o3_zcta_daily_parquet_manifest.csv`: row counts, date coverage, completeness, file sizes, and SHA-256 checksums.
- `o3_zcta_daily_qc_all_years.csv`: monthly completeness and distribution QC from the source aggregation.
- `README.md`: variables and citation notes.

Columns:

- `zip`
- `date`
- `pollutant`
- `year`
- `month`
- `temporal_resolution`
- `o3_ppb`
- `value_source`
- `fill_distance_m`
- `fill_cell`

The `o3_ppb` column is daily maximum 8-hour average ozone in parts per billion. `value_source`, `fill_distance_m`, and `fill_cell` provide audit metadata for nearest-cell fills; the completed 2005-2024 local outputs have zero missing ZCTA-day ozone values.

Citation:

Please cite the source ozone estimates using the Science article DOI supplied with these data: https://www.science.org/doi/10.1126/science.aed3197.

Also cite this repository/release for the ZCTA aggregation workflow. The reproducible aggregation script is `code/52_aggregate_daily_o3_to_zcta.R`, the QA mapping script is `code/53_map_daily_o3_zcta_examples.R`, and this packaging script is `code/54_package_daily_o3_release_assets.R`.
