# Air Pollution and Transplant Waitlist Outcomes

Reproduction code for the heart, kidney, liver, and lung waitlist study.
Primary models estimate cause-specific hazards of death or delisting for
deterioration, separately by organ and pollutant. The baseline models use
prelisting exposure. Time-varying models update exposure and organ-specific
severity at monthly intervals.

## Reproduce the paper

1. Use R 4.5 or later. Install the R packages listed in
   [the reproduction guide](docs/reproduction.md).
2. Obtain authorized access to SRTR SAFs. Set `SAF_Q1_2026_DIR` and
   `SAF_Q2_2025_DIR` as described in the guide.
3. From the repository root, run:

   ```sh
   Rscript code/run_paper.R --download
   ```

Run `Rscript code/run_paper.R --list` to see the exact order, or use
`Rscript code/run_paper.R displays validate` to regenerate figures and
tables from existing model outputs. A full run is computationally intensive.
Use a fresh checkout for a clean rebuild with no prior model checkpoints.

The generated [publication index](docs/reproduction.md#outputs) catalogs the
figures, tables, and supplementary analyses. Generated files, logs, and
candidate-level intermediates stay in the ignored `output/` directory.

## Repository layout

| Path | Contents |
| --- | --- |
| [code/run_paper.R](code/run_paper.R) | Ordered analysis and display workflow |
| [code/](code/) | Cohort construction, exposure linkage, and model fits |
| [code/figures_tables/](code/figures_tables/) | Figure and table generators |
| [code/revision_2026/](code/revision_2026/) | LAS, reviewer sensitivities, and publication postprocessing |
| [code/tests/](code/tests/) | Exposure-window and cohort logic checks |
| [data/reference/](data/reference/) | Public reference inputs and pinned asset checksums |
| [docs/reproduction.md](docs/reproduction.md) | Inputs, stages, and display index |

Environmental aggregation and release packaging are maintained in
[Climate Exposure Aggregation](https://github.com/petergraffy/climate-exposure-aggregation).
The separate [lung transplant survival-benefit study](https://github.com/petergraffy/lung_transplant_pollution_survival_benefit)
has its own repository.
