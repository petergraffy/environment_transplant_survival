# ZIP Crosswalk Plan, 2005-2024

This project uses the HUD-USPS ZIP Code Crosswalk files for year-aware ZIP-to-census geography assignment.

## Source

Primary source: HUD-USPS ZIP Code Crosswalk files.

HUD describes these files as quarterly snapshots derived from USPS ZIP+4 vacancy/address data. The ZIP-to-tract and ZIP-to-county files include residential, business, other, and total address ratios. For patient residential exposure assignment, `RES_RATIO` is the preferred weight.

Important limitations:

- HUD states that ZIP crosswalk files are not available before 2010 Q1.
- ZIP codes do not align cleanly with census, county, city, or state boundaries.
- PO Box-only ZIP codes are not included in the HUD files.
- 2010 Q1-2011 Q4 files use 2000 Census geographies.
- 2012 Q1-2022 Q4 files use 2010 Census geographies.
- 2023 Q1 onward files use 2020 Census geographies.

## Study-Year Strategy

- 2010-2024: use HUD Q4 files for each study year.
- 2005-2009: backcast the earliest HUD file, 2010 Q1, and flag these rows with `pre2010_backcast = TRUE`.

This is intentionally conservative and transparent. If the primary environmental exposure product begins in 2010 or later, prefer using the directly observed HUD years only. If the primary analysis begins in 2005, include a sensitivity analysis excluding 2005-2009.

## Outputs

The build script writes:

- `data/processed/crosswalks/zip_to_tract_crosswalk_2005_2024.parquet`
- `data/processed/crosswalks/zip_to_county_crosswalk_2005_2024.parquet`
- `data/processed/crosswalks/zip_crosswalk_2005_2024_qc.csv`
- `data/processed/crosswalks/zip_crosswalk_2005_2024_sources.csv`

Each crosswalk preserves multiple geography rows per ZIP per year and includes `res_ratio`, `bus_ratio`, `oth_ratio`, and `tot_ratio`.

## Raw Files

Place HUD files in:

`data/raw/hud_usps_crosswalk/`

Expected naming pattern:

- `ZIP_TRACT_032010.xlsx`
- `ZIP_COUNTY_032010.xlsx`
- `ZIP_TRACT_122010.xlsx`
- `ZIP_COUNTY_122010.xlsx`
- ...
- `ZIP_TRACT_122024.xlsx`
- `ZIP_COUNTY_122024.xlsx`

The script also accepts `.xls`, `.csv`, or `.txt` with the same stem.

Run:

```sh
Rscript code/exploratory/01_build_zip_crosswalk_2005_2024.R
```

The script has an optional download mode:

```sh
Rscript code/exploratory/01_build_zip_crosswalk_2005_2024.R --download
```

HUD may block command-line downloads with a web challenge. If so, download the files manually through HUD USER and rerun without `--download`.
