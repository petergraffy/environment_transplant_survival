# Reproduction guide

## Inputs

- **Restricted SRTR data.** Set `SAF_Q1_2026_DIR` to the directory containing
  `pubsaf2603/` and `Supplement2603/`. The official LAS audit uses
  `Lung_supplement2603/las_saf2603.sas7bdat` and `cas_2603.sas7bdat` from
  that directory. Set `SAF_LUNG_SUPPLEMENT_DIR` for a separate location.
  Set `SAF_Q2_2025_DIR` for the historical center-map reference, which uses
  the Q2 2025 SAF by default. Use only data you are authorized to access.
- **Public exposure data and boundaries.** `data/reference/publication_inputs.csv`
  lists the 50 Parquet assets and the 2020 Census ZCTA archive, with URLs and
  SHA-256 hashes. `Rscript code/run_paper.R inputs --download` retrieves missing
  files, verifies them, and extracts the boundaries. A checksum mismatch stops
  the build so an old or substituted input cannot silently enter the analysis.
- **Reference tables.** The ACS covariates, 2022 CDC SVI ZCTA data, and US-CRS
  point mapping are included. The optional `code/06_build_community_covariates.R`
  documents the ACS extraction and requires `CENSUS_API_KEY` to regenerate it.

See [exposure definitions and coverage](daily_pollution_analysis_inputs.md).
The pinned annual NO2 asset includes the latitude correction for 2025. The
versioned monthly release starts January 2019.

## Software

Use R 4.5 or later with `arrow`, `broom`, `data.table`, `dplyr`, `fixest`,
`ggplot2`, `gt`, `haven`, `maps`, `patchwork`, `readr`, `scales`, `sf`, `stringr`,
`survival`, `tibble`, `tidyr`, `tidyselect`, and `viridis`. To rebuild the ACS
reference, also install `jsonlite`. The optional PDF figure QA stage uses Python with
`pypdfium2` and `Pillow`. Set `PAPER_PYTHON` for a custom Python path.

## Run order

From the repository root, `Rscript code/run_paper.R --list` prints the
script order, and `--check` checks that all entry points exist and that R
sources parse. With inputs already on disk, run `Rscript code/run_paper.R`.
Run named stages individually when upstream inputs already exist.

| Stage | Work |
| --- | --- |
| `inputs` | Verify exposure assets and boundaries. Add `--download` to fetch missing assets. |
| `cohort` | Prepare the deduplicated cohort, daily/monthly pollution caches, registration intervals, and center geography. |
| `models` | Fit the reference baseline and monthly time-varying models, severity associations, and subgroup and multipollutant analyses. Compute reference AJ curves and quartiles. |
| `sensitivities` | Fit and audit PH, listing eras, registration spells, concurrency, death-only outcomes, and the matched 2015-2019 Fine-Gray comparison. Validate the ACS vulnerability proxy. |
| `publication` | Fit the official LAS models, update the lung rows in the shared results, and run the extended lung follow-up and transplant/improvement analyses. |
| `displays` | Render the Cox forests, AJ curves and risk tables, flow diagram, geography and risk maps, Table 1, and all supplemental tables and reports. |
| `validate` | Reconcile publication model estimates and denominators. Export HTML tables and build the display index and manifest. |
| `qa` | Optional. Render publication PDFs and create visual review contact sheets. |

The publication models restrict lung entry and follow-up through March 8,
2023 and use official LAS at monthly interval starts. The later cohort-flow
and transplant/improvement scripts instead use the separately saved extended
follow-up analysis. That analysis caps LAS updates at the transition and
carries the last observed LAS forward. The runner preserves both definitions.

Use a fresh checkout for a full rebuild. Several legacy analysis steps reuse
private caches or per-model checkpoints under `output/`. When the cohort or
inputs change, a clean checkout prevents stale results from entering a run.
The full build is sequential and can take hours on the national cohort.

## Outputs

- `output/revision_2026_checklist/las_era_lung/publication_index.html` links
  the figures, tables, and supplements. `publication_manifest.csv` lists the
  same artifacts, and `publication_validation.txt` records final checks.
- `output/figures/` and `output/tables/` hold the manuscript displays and
  their source summaries. Reviewer diagnostics and reports live in
  `output/revision_2026_checklist/`.
- `output/paper_run_logs/` contains one log per script, grouped by invocation.
  `Rscript code/run_paper.R displays validate` rebuilds the displays and index
  from existing model results. Use `Rscript code/run_paper.R qa` for PDF QA.

Generated files and candidate-level analysis caches are private working files
and are excluded from Git.
