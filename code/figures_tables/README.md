# Figures and tables

Run `Rscript code/run_paper.R displays validate` from the repository root
after the model and publication stages. The validation stage writes
`output/revision_2026_checklist/las_era_lung/publication_index.html` and a
CSV manifest. Figure scripts emit their analysis-specific PDF/PNG and source
summaries. The runner lists their execution order with `--list`.

- 19 prepares the center location tables used by later maps.
- 60-62 and 92-94 create concentration, cohort geography, and modeled-risk maps.
- 69, 75, 77, and 89 create descriptive and sensitivity tables.
- 71 and 83 create Aalen-Johansen displays.
- 76 invokes the cohort-flow generator, which the runner also lists directly.
- 82, 85-87, and 91 create model forest plots.
- The [reviewer analysis and postprocessing scripts](../revision_2026/) create
  risk tables, diagnostics, the missingness table, and supplement reports.
