# Daily LGHAP PM2.5 by CONUS county

The county pipeline aggregates the available 2005-2024 daily LGHAP PM2.5
NetCDF archives directly to fixed 2020 Census county boundaries. It includes
3,108 counties and county equivalents in the contiguous 48 states and DC.
Alaska, Hawaii, and territories are excluded. Earlier source years are not
present in the local archive inventory.

## Method

`code/pollution_aggregation/69_aggregate_lghap_daily_pm25_to_county.R` uses
`exactextractr` to calculate each grid cell's county overlap fraction and
cell area. Their product is the spatial weight. A cached sparse matrix
applies these weights to each daily PM2.5 field. The daily estimate is
the weighted sum of valid concentrations divided by their total weight.
This is an area-weighted ambient concentration, not a population-weighted
exposure or a mean of previously aggregated ZIP values.

`valid_area_fraction` records the proportion of the county's mapped area
with valid source values that day. Missing source cells are excluded from
the weighted denominator, not assigned zero. A county with no valid source
coverage is flagged and filled from the nearest valid raster cell within
500 km of a point on its surface, matching the existing ZCTA fallback.
`value_source = nearest_raster_cell`, `fill_distance_m`, and `fill_cell`
identify these exceptions, while `valid_area_fraction` remains zero.
An unfillable county causes the job to stop. NetCDF fill values and invalid
negative values are checked.

The grid's coordinate vectors are checked against the cached grid every day.
On the first day of every month, the indexed sparse calculation is checked
against independent raster-based area-weighted extraction for Los Angeles
County (06037), Cook County (17031), and New York County (36061).

County boundaries are the Census Bureau's
[2020 1:500,000 cartographic boundary file](https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_us_county_500k.zip).
All years use this same geography, including Connecticut's 2020 counties.
These simplified cartographic boundaries are not full-resolution TIGER/Line
boundaries. They match the boundary vintage and scale of the existing ZCTA
workflow, although this county pipeline uses exact area weights rather than
the ZCTA pipeline's rasterized zonal means.

## Outputs and joins

The GitHub release `lghap-pm25-county-conus-v1` contains one daily Parquet
file per year (`lghap_pm25_county_daily_YYYY.parquet`) and one annual Parquet
file (`lghap_pm25_county_annual.parquet`). Both store FIPS codes as strings.
`lghap_pm25_county_parquet_manifest.csv` records coverage, row counts, and
SHA-256 checksums for the Parquet files. The release also includes the
monthly QC, county lookup, and fill audit listed below.

Run `Rscript code/pollution_aggregation/70_package_lghap_daily_pm25_county_release_assets.R`
from the repository root to build and verify the release assets in
`data/release/lghap_pm25_county_parquet/`.

Outputs are in `data/processed/lghap_pm25_county_daily/`:

- `lghap_pm25_county_daily_YYYY_MM.csv.gz`: one row per county and day.
- `lghap_pm25_county_daily_YYYY_MM_qc.csv`: monthly QC and SHA-256 checksum.
- `lghap_pm25_county_daily_manifest.csv`: combined monthly QC.
- `lghap_pm25_county_annual.csv`: arithmetic mean of all daily values in each year; leap days included.
- `county_lookup_2020_conus.csv`: county FIPS codes and names.
- `lghap_pm25_county_fill_audit.csv`: every nearest-cell fill and its distance.

Join daily data using `county_fips` and `date`, or annual data using
`county_fips` and `year`. Read FIPS codes as strings to preserve leading zeros.
`pm25_ug_m3` is PM2.5 in micrograms per cubic meter. Daily outputs also contain
`state_fips`, `county_name`, `year`, `month`, `pollutant`,
`temporal_resolution`, `value_source`, `valid_area_fraction`,
`fill_distance_m`, and `fill_cell`. Annual outputs retain `nearest_fill_days`.

## Running and progress

### Completed-run validation

The 2005-2024 run completed on November 20, 2026: 240 monthly files,
7,305 dates, 22,703,940 county-day rows, and 62,160 county-year rows.
All counties have every expected date, with no duplicate county/date keys
or missing final concentrations. All monthly SHA-256 checksums were verified.

One source gap required a nearest-cell fill: Richland County, Wisconsin
(55103), on May 2, 2023. There were no valid source grid cells inside the
county on that date. The filled value is 3 micrograms per cubic meter from
a valid grid cell 47,082.49 m from the county's interior point. It is retained
in the fill audit and counted in the county's 2023 annual `nearest_fill_days`.
All other county-day values are direct area-weighted source-grid averages.

### Command

```powershell
Rscript code/pollution_aggregation/69_aggregate_lghap_daily_pm25_to_county.R
```

Optional environment variables: `LGHAP_COUNTY_START_YEAR` (2005),
`LGHAP_COUNTY_END_YEAR` (2024), `LGHAP_COUNTY_WORKERS` (12), and
`LGHAP_COUNTY_DOWNLOADS` (the user's Downloads directory).
The repository's `data/raw/lghap_pm25_daily` is searched first.
The 2020 county shapefile must be present at
`data/cache/cb_2020_us_county_500k/cb_2020_us_county_500k.shp`.

Years run sequentially, with monthly tasks assigned to the worker pool.
The worker log is `logs/lghap_pm25_county_workers.log`. Monthly outputs are
checkpointed and can be resumed after checking their checksums and QC.
Each worker extracts one daily file at a time to its own temporary directory
and removes the extracted file after aggregation. Source archives are retained.
The final step rereads all outputs and checks dates, duplicate keys, finite
concentrations, and complete daily coverage for every county.

## Citation

Bai, K., Li, K., Shao, L., Li, X., Liu, C., Li, Z., Ma, M., Han, D., Sun, Y.,
Zheng, Z., Li, R., Chang, N.-B., and Guo, J. (2024). LGHAP v2: a global gap-free
aerosol optical depth and PM2.5 concentration dataset since 2000 derived via
big Earth data analytics. *Earth System Science Data*, 16, 2425-2448.
https://doi.org/10.5194/essd-16-2425-2024

The 2022-2024 inputs are downloaded LGHAP PM2.5 products, not the separate
locally derived AOD-to-PM2.5 estimates. Also cite the Census boundaries and
this repository's county aggregation workflow when using these outputs.
