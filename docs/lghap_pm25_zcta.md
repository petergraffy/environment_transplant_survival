# LGHAP Global PM2.5 ZCTA Aggregates

`code/pollution_aggregation/39_aggregate_lghap_pm25_to_zcta.R` aggregates the downloaded LGHAP
Global PM2.5 NetCDF archive to 2020 ZCTA5 polygons.

The downloaded archive inspected in `C:/Users/Peter Graffy/Downloads/11435542.zip`
contains one file per year from 2000-2021:

`LGHAP.Global_PM25.Y001.AYYYY.nc`

Each NetCDF has one `PM25` layer with dimensions `lat x lon`; there is no
`time` or date dimension. The derived ZCTA output is therefore annual, not
daily.

## Outputs

Yearly ZCTA files:

`data/processed/lghap_pm25_zcta_annual/lghap_pm25_zcta_annual_YYYY.csv.gz`

Manifest and QA:

- `data/processed/lghap_pm25_zcta_annual/lghap_pm25_zcta_annual_manifest.csv`
- `data/processed/lghap_pm25_zcta_annual/lghap_pm25_zcta_annual_qc_all_years.csv`

## Columns

- `zip`: 2020 ZCTA5 code
- `pollutant`: `pm25_lghap`
- `year`: calendar year
- `temporal_resolution`: `annual`
- `pm25_ug_m3`: ZCTA mean PM2.5 concentration
- `value_source`: `zonal_mean`, `nearest_raster_cell`, or `missing`
- `fill_distance_m`: distance for nearest-cell fills
- `fill_cell`: source raster cell for nearest-cell fills

ZCTAs without a direct zonal grid value are filled from the nearest valid raster
cell, matching the audit pattern used in the other ZCTA exposure pipelines.
