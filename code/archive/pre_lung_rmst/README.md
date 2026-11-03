# Historical code before the primary lung RMST analysis

Archived November 2, 2026, when the lung transplant RMST workflow became the
primary analysis. This archive preserves the previous four-organ waitlist
models, publication figures, revision analyses, exploratory work and their tests.
Statements about the "primary analysis" inside these historical files describe
the earlier study, not the current specification.

The folder layout under the former `code/` directory is retained here:

- Root numbered scripts: earlier waitlist models and sensitivity analyses.
- `revision_2026/`: manuscript revision and official LAS-era preparation.
- `figures_tables/`: earlier publication graphics and tables.
- `exploratory/`: previous feasibility, diagnostic and exploratory analyses.
- `tests/`: tests specific to those historical analyses.

Run archived scripts from the **repository root**. Source/script references have
been changed to `code/archive/pre_lung_rmst/...`; shared input helpers remain
under `code/`. Data and historical outputs retain their original paths.

The current cohort-preparation adapter still reuses the established cohort
functions and `revision_2026/30_prepare_las_era_lung.R`. Archiving does not imply
these provenance dependencies can be deleted. It separates historical analysis
entry points from the active lung RMST workflow without changing cohort logic.

The tracked-file relocation inventory is in `manifest.csv`. Earlier committed
content is preserved by Git history. Local pre-existing, uncommitted edits are
retained at their relocated working-tree paths and are not bundled into the
primary-analysis commit.

Start current work with [`code/run_primary_lung_analysis.R`](../../run_primary_lung_analysis.R)
and the [primary methods](../../../docs/lung_pollution_survival_benefit.md).
