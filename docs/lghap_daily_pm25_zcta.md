# LGHAP Daily PM2.5 ZCTA Aggregates

`code/pollution_aggregation/40_aggregate_lghap_daily_pm25_to_zcta.R` aggregates the downloaded LGHAP
daily PM2.5 archive to 2020 ZCTA5 polygons for CONUS.

The 2021 archive inspected at `C:/Users/Peter Graffy/Downloads/8313615.zip`
contains one nested monthly zip per month:

`LGHAP.Global_PM25.D001.M2021MM.zip`

Each monthly zip contains one NetCDF per day:

`LGHAP.Global_PM25.D001.AYYYYMMDD.nc`

Each daily NetCDF has one `PM25` layer with dimensions `lat x lon`; the script
corrects that orientation before aggregating to ZCTA.

## Outputs

Daily ZCTA rows are written one month per file:

`data/processed/lghap_pm25_zcta_daily/lghap_pm25_zcta_daily_2021_MM.csv.gz`

Manifest and QA:

- `data/processed/lghap_pm25_zcta_daily/lghap_pm25_zcta_daily_manifest_2021.csv`
- `data/processed/lghap_pm25_zcta_daily/lghap_pm25_zcta_daily_qc_2021_all_months.csv`

## Columns

- `zip`: 2020 ZCTA5 code
- `date`: daily date
- `pollutant`: `pm25_lghap_daily`
- `year`: calendar year
- `month`: calendar month
- `temporal_resolution`: `daily`
- `pm25_ug_m3`: daily ZCTA mean PM2.5 concentration
- `value_source`: `zonal_mean`, `nearest_raster_cell`, or `missing`
- `fill_distance_m`: distance for nearest-cell fills
- `fill_cell`: source raster cell for nearest-cell fills

ZCTAs without a direct zonal grid value are filled from the nearest valid raster
cell, matching the audit pattern used in the other ZCTA exposure pipelines.
