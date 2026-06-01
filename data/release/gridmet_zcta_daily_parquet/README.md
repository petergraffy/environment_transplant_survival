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

## Source Data Citation

Please cite the original gridMET data source when using these derived ZCTA
aggregates:

- Abatzoglou, J. T. (2013). Development of gridded surface meteorological data
  for ecological applications and modelling. *International Journal of
  Climatology*, 33(1), 121-131. https://doi.org/10.1002/joc.3413
- gridMET data portal: https://www.northwestknowledge.net/metdata/data/
