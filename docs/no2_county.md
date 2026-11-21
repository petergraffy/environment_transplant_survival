Annual and monthly surface NO2 for 3,108 counties and county equivalents in the contiguous U.S. and DC. Concentrations are area-weighted over fixed 2020 Census county boundaries and reported in `no2_ppbv` (parts per billion by volume).

- `no2_county_annual.parquet` covers 2005-2025. Join on `county_fips` and `year`. The annual source switches from Mohegh and Anenberg after 2020 to Nawaz TROPOMI v2 in 2021.
- `no2_county_monthly.parquet` covers May 2018 through December 2025 using Nawaz TROPOMI v2. Join on `county_fips`, `year`, and `month`.

`county_fips` is a five-character string. `source_product` and `source_file` identify each input. Annual TROPOMI grid rows are mapped using the NetCDF latitude coordinates. `valid_area_fraction` reports valid source coverage. The QC and fill audit identify any nearest-cell fills.

**Citations**

- Mohegh, A., and Anenberg, S. (2021). *Global surface NO2 concentrations 1990-2020*. [Figshare dataset](https://doi.org/10.6084/m9.figshare.12968114) and [Anenberg et al. (2022)](https://doi.org/10.1016/S2542-5196(21)00255-2).
- Nawaz, M. O. (2026). *Monthly and Annual US TROPOMI Surface NO2 Estimates (~1 km x 1 km)*, v2. [Zenodo dataset](https://doi.org/10.5281/zenodo.18919769) and [Nawaz et al. (2025)](https://doi.org/10.1021/acsestair.4c00153).
