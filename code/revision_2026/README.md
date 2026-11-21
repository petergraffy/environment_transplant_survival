# Publication analyses and postprocessing

The [paper runner](../run_paper.R) is the ordered entry point. It prepares
registration-level sensitivities and diagnostics, builds official LAS histories,
fits the lung publication models, and combines those estimates with the other
organs before rendering the displays.

- 29-31 read official LAS, prepare the lung cohort, and fit its models.
- 32 synchronizes the lung results and AJ curves. It uses the original
  reference-cohort quartile thresholds and verifies the other organs' rows.
- 33 exports outcome counts and the model covariate missingness table.
- 35 and 37 update diagnostics and severity associations for official LAS.
- 40-43 retain the extended lung follow-up report, cohort flow, and baseline
  transplant/improvement sensitivity. The extended bundle keeps the last LAS
  observed by March 8, 2023 and has its own result files.
- 39 verifies denominators and estimates, exports copyable HTML tables, and
  creates the publication index and manifest. 38 renders PDF previews for QA.

The publication model bundle censors lung follow-up at the LAS-to-CAS
transition. The later cohort-flow and transplant/improvement analyses use the
extended-follow-up bundle, as in the existing figure and supplement scripts.
See [the reproduction guide](../../docs/reproduction.md) for inputs and stages.
