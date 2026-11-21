Daily and annual surface ozone for 3,108 counties and county equivalents in the contiguous U.S. and DC, 2005-2024. Daily maximum 8-hour average (MDA8) ozone is area-weighted over fixed 2020 Census county boundaries. Annual values average all daily values in each year.

- `o3_county_daily_YYYY.parquet` has one row per county and day. Join on `county_fips` and `date`.
- `o3_county_annual.parquet` has one row per county and year. Join on `county_fips` and `year`.

`county_fips` is a five-character string. `o3_ppb` is in parts per billion. The release includes the county lookup, monthly QC, Parquet SHA-256 checksums, and a nearest-cell fill audit. `valid_area_fraction` gives the share of each county with valid source data. Any zero-coverage day filled from a nearby grid cell is flagged in `value_source`.

**Citation**

Deng, W. et al. (2026). *Fires reverse progress toward ozone air quality standards in the United States*. Science, 392, 1088-1092. [Article](https://doi.org/10.1126/science.aed3197) | [Source dataset](https://doi.org/10.5281/zenodo.17502917)
