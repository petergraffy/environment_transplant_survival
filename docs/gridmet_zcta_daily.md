# gridMET Daily ZCTA Weather Exposures

`code/20_aggregate_gridmet_to_zcta.R` downloads daily gridMET NetCDF files from
`https://www.northwestknowledge.net/metdata/data/` and aggregates them to 2020
ZCTA5 polygons.

## Outputs

Daily wide files are split by calendar year:

`data/processed/gridmet_zcta_daily/gridmet_zcta_daily_YYYY.csv.gz`

Combined metadata:

- `data/processed/gridmet_zcta_daily/gridmet_zcta_daily_manifest_all_years.csv`
- `data/processed/gridmet_zcta_daily/gridmet_zcta_daily_qc_all_years.csv`

Raw NetCDF downloads are cached in:

`data/raw/gridmet/`

## Variables

- `tmax_c`: daily maximum temperature, degrees C
- `tmin_c`: daily minimum temperature, degrees C
- `tmean_c`: mean of `tmax_c` and `tmin_c`, degrees C
- `temp_range_c`: `tmax_c - tmin_c`, degrees C
- `rmax_pct`: daily maximum relative humidity, percent
- `rmin_pct`: daily minimum relative humidity, percent
- `rhmean_pct`: mean of `rmax_pct` and `rmin_pct`, percent
- `pr`: daily precipitation, mm/day
- `sph`: specific humidity, kg/kg
- `srad`: downward shortwave radiation
- `vs`: wind speed
- `pet`: potential evapotranspiration
- `etr`: reference evapotranspiration
- `erc`: energy release component fire-danger index

Each source variable also has `_source`, `_fill_distance_m`, and `_fill_cell`
audit columns. ZCTAs without a direct zonal grid value are filled from the
nearest valid gridMET raster cell.
