param([string]$RunDirectory = 'output/revision_runs/rolling_20260922')
$ErrorActionPreference = 'Stop'
$root = (Get-Location).Path
$run = Join-Path $root $RunDirectory
$before = Join-Path $run 'before'
New-Item -ItemType Directory -Force -Path $run | Out-Null
if (-not (Test-Path -LiteralPath (Join-Path $run 'snapshot_complete.txt'))) {
    $folders = @(
        'prior_year_pollution_cox_svi', 'prior_year_initial_listing_cox_model_set',
        'prior_year_pollution_baseline_severity_associations', 'prior_year_pollution_quartile_aalen_johansen_cif',
        'timevarying_pollution_cox_svi', 'timevarying_pollution_severity_models_with_listing_year',
        'prior_year_and_timevarying_pollution_subgroup_cox', 'updated_methods_sensitivity_models',
        'tables', 'figures/figure1_cox_baseline_timevarying',
        'figures/prior_year_initial_listing_cox_model_set', 'figures/pollution_severity_associations',
        'figures/pollution_severity_associations_with_listing_year',
        'figures/prior_year_pollution_quartile_aalen_johansen_cif',
        'figures/prior_year_pollution_quartile_aalen_johansen_simplified',
        'figures/prior_year_and_timevarying_pollution_subgroup_cox',
        'figures/pollution_absolute_1yr_survival_maps', 'figures/candidates_heart_lung_absolute_risk_abc',
        'figures/primary_pollution_study_period_maps', 'figures/transplant_centers_by_organ_four_panel'
    )
    foreach ($folder in $folders) {
        $source = Join-Path $root "output/$folder"
        if (-not (Test-Path -LiteralPath $source)) { continue }
        foreach ($file in Get-ChildItem -LiteralPath $source -File -Recurse) {
            if ($file.FullName -match '[\\/]cache[\\/]') { continue }
            $relative = [System.IO.Path]::GetRelativePath((Join-Path $root 'output'), $file.FullName)
            $target = Join-Path $before $relative
            New-Item -ItemType Directory -Force -Path (Split-Path $target) | Out-Null
            Copy-Item -LiteralPath $file.FullName -Destination $target
        }
    }
    Copy-Item -LiteralPath (Join-Path $root 'code') -Destination (Join-Path $run 'code_at_start') -Recurse
    Set-Content -LiteralPath (Join-Path $run 'snapshot_complete.txt') -Value (Get-Date -Format o)
}
$steps = @(
    'code/81_prior_year_pollution_cox_svi.R',
    'code/90_prior_year_initial_listing_cox_model_set.R',
    'code/80_prior_year_pollution_baseline_severity_associations.R',
    'code/78_prior_year_pollution_quartile_aalen_johansen_cif.R',
    'code/79_timevarying_pollution_cox_svi.R',
    'code/84_timevarying_pollution_severity_models.R',
    'code/86_prior_year_and_timevarying_pollution_subgroup_cox.R',
    'code/88_updated_methods_sensitivity_models.R',
    'code/figures_tables/82_plot_figure1_prior_year_and_timevarying_cox.R',
    'code/figures_tables/83_plot_prior_year_aalen_johansen_simplified.R',
    'code/figures_tables/85_plot_pollution_severity_associations.R',
    'code/figures_tables/86_plot_prior_year_and_timevarying_pollution_subgroup_forest.R',
    'code/figures_tables/87_plot_heart_lung_pm25_no2_subgroup_forest.R',
    'code/figures_tables/89_make_updated_sensitivity_wide_tables.R',
    'code/figures_tables/91_plot_initial_listing_model_tier_forest.R',
    'code/figures_tables/60_plot_primary_pollution_study_period_maps.R',
    'code/figures_tables/69_make_table1_baseline_characteristics.R',
    'code/figures_tables/92_plot_pollution_absolute_1yr_survival_maps.R',
    'code/figures_tables/93_plot_candidates_heart_lung_absolute_risk_abc.R',
    'code/figures_tables/94_plot_transplant_centers_by_organ_four_panel.R'
)
foreach ($script in $steps) {
    $name = [System.IO.Path]::GetFileNameWithoutExtension($script)
    $done = Join-Path $run "$name.done"
    if (Test-Path -LiteralPath $done) { continue }
    $log = Join-Path $run "$name.log"
    Write-Output "$(Get-Date -Format o) RUNNING $script"
    & Rscript $script *> $log
    $exitCode = $LASTEXITCODE
    if ($exitCode -ne 0) {
        Get-Content -LiteralPath $log -Tail 30
        throw "Step failed ($exitCode): $script. Completed steps can be resumed."
    }
    Set-Content -LiteralPath $done -Value (Get-Date -Format o)
    Write-Output "$(Get-Date -Format o) COMPLETE $script"
}
Write-Output 'All revision model and figure steps completed.'
