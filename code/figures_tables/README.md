# Figure and Table Scripts

This folder contains scripts that create manuscript figures, supplemental
figures, and tables from model outputs.

Run these scripts from the repository root. Most scripts read from `output/`
analysis products and write to `output/figures/` or `output/tables/`.

Key manuscript assets:

- `59_plot_primary_waitlist_pollution_cox_figures.R`: primary Cox forest plots.
- `60_plot_primary_pollution_study_period_maps.R`: study-period pollution maps.
- `61_plot_waitlist_candidates_zcta_map.R` and
  `62_plot_waitlist_geography_ab_maps.R`: candidate and transplant-center maps.
- `69_make_table1_baseline_characteristics.R`: Table 1.
- `70_plot_primary_pollution_quartile_aalen_johansen_cif.R`: primary
  Aalen-Johansen cumulative incidence curves.
- `71_plot_heart_lvad_kidney_dialysis_aalen_johansen_cif.R`: subgroup
  Aalen-Johansen curves for heart durable LVAD and kidney dialysis duration.
- `73_make_sensitivity_supplement_tables.R`: supplemental sensitivity tables.
- `76_make_cohort_flow_diagram.R`: cohort flow diagram.
