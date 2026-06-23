# Daily Ozone ZCTA Aggregates

`code/pollution_aggregation/52_aggregate_daily_o3_to_zcta.R` aggregates downloaded daily CONUS ozone
NetCDF archives to 2020 ZCTA5 polygons.

The inspected 2024 archive at `C:/Users/Peter Graffy/Downloads/2024.zip`
contains one NetCDF per day:

`YYYY/YYYYMMDD.nc`

Each daily NetCDF has one `mda8` layer in ppb. The script writes one compressed
CSV per year-month:

`data/processed/o3_zcta_daily/o3_zcta_daily_YYYY_MM.csv.gz`

## Columns

- `zip`: 2020 ZCTA5 code
- `date`: daily date
- `pollutant`: `o3_daily_mda8`
- `year`: calendar year
- `month`: calendar month
- `temporal_resolution`: `daily`
- `o3_ppb`: daily MDA8 ozone in ppb
- `value_source`: `zonal_mean`, `nearest_raster_cell`, or `missing`
- `fill_distance_m`: distance for nearest-cell fills
- `fill_cell`: source raster cell for nearest-cell fills

## Example

```powershell
Set-Item -Path Env:DAILY_O3_YEAR -Value '2024'
Set-Item -Path Env:DAILY_O3_ZIP -Value 'C:/Users/Peter Graffy/Downloads/2024.zip'
Rscript code\52_aggregate_daily_o3_to_zcta.R
```

The output includes monthly QC files and an all-month yearly QC summary.

## Release Packaging

`code/pollution_aggregation/54_package_daily_o3_release_assets.R` converts monthly CSV outputs to
yearly Parquet files in `data/release/o3_zcta_daily_parquet`. The release
assets include yearly Parquet files, a manifest with SHA-256 checksums, all-year
QC, and a README.

`code/pollution_aggregation/55_upload_daily_o3_release_assets.ps1` creates or updates the GitHub
release tagged `o3-zcta-daily-v1`.

## Citation

Please cite the source ozone estimates using the Science article DOI supplied
with these data:

https://www.science.org/doi/10.1126/science.aed3197
