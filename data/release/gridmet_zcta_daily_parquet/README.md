# gridMET ZCTA Daily Parquet Release Assets

These files contain value-only daily gridMET exposures aggregated to 2020 ZCTA5 polygons.

Each yearly Parquet file contains:

- `zip`
- `date`
- `year`
- `tmax_c`
- `tmin_c`
- `tmean_c`
- `temp_range_c`
- `rmax_pct`
- `rmin_pct`
- `rhmean_pct`
- `pr`
- `sph`
- `srad`
- `vs`
- `pet`
- `etr`
- `erc`

Audit/fill summaries are in `gridmet_zcta_daily_qc_all_years.csv` and the committed documentation.
The original local processing outputs include source/fill columns and can be regenerated with `code/20_aggregate_gridmet_to_zcta.R`.
