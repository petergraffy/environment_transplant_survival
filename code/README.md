# Analysis code

Run from the repository root. [run_paper.R](run_paper.R) is the execution
order for the cohort, models, sensitivities, figures, tables, and validation.
The numbered filenames preserve links to the publication-era scripts.

- `50_primary_waitlist_period_pollution_cox.R` prepares the deduplicated input
  cohort. Its historical filename and output path are retained for compatibility.
- Scripts 78-90 produce the rolling-prelisting and monthly time-varying results.
- [revision_2026/](revision_2026/) contains LAS preparation, reviewer analyses,
  publication updates, and report export.
- [figures_tables/](figures_tables/) contains display generators.
- [tests/](tests/) contains reproducibility checks.

See [inputs, stages, and outputs](../docs/reproduction.md).
