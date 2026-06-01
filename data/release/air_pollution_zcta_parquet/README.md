# Air Pollution ZCTA Parquet Release Assets

These files contain final filled air-pollution exposure surfaces aggregated to 2020 ZCTA5 polygons.

Assets:

- `air_pollution_zcta_pm25_monthly_2005_2023.parquet`: monthly PM2.5, `pm25_ug_m3`.
- `air_pollution_zcta_o3_monthly_2005_2023.parquet`: monthly ozone, `o3_ppb`.
- `air_pollution_zcta_no2_annual_2005_2024.parquet`: annual NO2, `no2`.
- `air_pollution_zcta_parquet_manifest.csv`: sizes, row counts, completeness, and SHA-256 checksums.
- `fill_summary.csv`: nearest-fill audit by source file.
- `fill_audit_summary.csv`: additional fill audit summary when available.

Columns include `zip`, `pollutant`, `year`, `month`, `temporal_resolution`, the pollutant value column, `value_source`, and `fill_distance_km`.
PM2.5 and ozone currently cover 2005-2023 in the local source files; NO2 covers 2005-2024.

## Source Data Citations

Please cite the original pollutant data sources when using these derived ZCTA
aggregates:

- PM2.5: Shen, S., Li, C., van Donkelaar, A., Jacobs, N., Wang, C., & Martin,
  R. V. (2024). Enhancing global estimation of fine particulate matter
  concentrations by including geophysical a priori information in deep
  learning. *ACS ES&T Air*. https://doi.org/10.1021/acsestair.3c00054
- PM2.5 data portal: WashU Atmospheric Composition Analysis Group SatPM2.5,
  https://sites.wustl.edu/acag/surface-pm2-5/
- Ozone: Liu, R., Chu, L., Deziel, N. C., & Chen, K. (2026). Four-decade
  (1980-2023) surface ozone concentrations across the contiguous United States:
  Fine-resolution estimates and health implications. *Environmental Science &
  Technology*. https://doi.org/10.1021/acs.est.5c16412
- Ozone data DOI: Yale Dataverse, https://doi.org/10.60600/YU/M1WT9R
- NO2: Mohegh, A., & Anenberg, S. (2021). Global surface NO2 concentrations
  1990-2020. Figshare. https://doi.org/10.6084/m9.figshare.12968114
- NO2: Nawaz, M. O. (2025). Monthly and annual US TROPOMI surface NO2 estimates
  (~1 km x 1 km), version 1.01. Zenodo.
  https://doi.org/10.5281/zenodo.14646034
