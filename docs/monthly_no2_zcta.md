# Monthly NO2 ZCTA Aggregates

`code/pollution_aggregation/56_aggregate_monthly_no2_to_zcta.R` aggregates monthly TROPOMI/LUR
CONUS surface NO2 NetCDF4 files to 2020 ZCTA5 polygons.

The local archive used for the 2019-2025 run was:

`C:/Users/Peter Graffy/Downloads/18919769.zip`

It contains monthly files named like:

`monthly_mean_tropomi_lur_conus_surface_no2_MMYYYY.v1.02_03082026.nc4`

Each NetCDF4 file has one `Surface_NO2` layer in `ppbv`. The script writes one
compressed CSV per year-month:

`data/processed/no2_zcta_monthly/no2_zcta_monthly_YYYY_MM.csv.gz`

## Columns

- `zip`: 2020 ZCTA5 code
- `pollutant`: `no2`
- `year`: calendar year
- `month`: calendar month
- `temporal_resolution`: `monthly`
- `no2_ppbv`: monthly surface NO2 in ppbv
- `value_source`: `zonal_mean`, `nearest_raster_cell`, or `missing`
- `fill_distance_m`: distance for nearest-cell fills
- `fill_cell`: source raster cell for nearest-cell fills

## Example

```powershell
Set-Item -Path Env:MONTHLY_NO2_YEARS -Value '2019:2025'
Set-Item -Path Env:MONTHLY_NO2_ZIP -Value 'C:/Users/Peter Graffy/Downloads/18919769.zip'
Set-Item -Path Env:MONTHLY_NO2_SKIP_EXISTING_MONTH -Value 'true'
Rscript code\56_aggregate_monthly_no2_to_zcta.R
```

The output includes monthly QC files and an all-month QC summary.

## Release Packaging

`code/pollution_aggregation/57_package_monthly_no2_release_assets.R` converts monthly CSV outputs to
yearly Parquet files in `data/release/no2_zcta_monthly_parquet`. The release
assets include yearly Parquet files, a manifest with SHA-256 checksums, all-year
QC, and a README.

`code/pollution_aggregation/58_upload_monthly_no2_release_assets.ps1` creates or updates the GitHub
release tagged `no2-zcta-monthly-v1`.

## Citation

Nawaz, M. Omar. Monthly and Annual US TROPOMI Surface NO2 Estimates (~1km x
1km), Version v2. Zenodo, 2026. https://doi.org/10.5281/zenodo.18919769.

Zenodo record: https://zenodo.org/records/18919769.

Related publication DOI listed by Zenodo:
https://doi.org/10.1021/acsestair.4c00153.

Also cite this repository/release for the ZCTA aggregation workflow.
