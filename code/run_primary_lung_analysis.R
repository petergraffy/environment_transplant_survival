#!/usr/bin/env Rscript
# Run each stage in a fresh R process so model data do not accumulate in memory.
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% c("--reuse-fitted", "--help"))) stop("Unknown option. Use --help.")
if ("--help" %in% args) {
  cat("Primary lung RMST analysis\n",
    "  Rscript code/run_primary_lung_analysis.R\n",
    "  Rscript code/run_primary_lung_analysis.R --reuse-fitted\n",
    "Full run: verify/prepare cohort, refit models, calculate 365-day RMST, bootstrap and validate.\n",
    "--reuse-fitted: regenerate RMST and restore intervals using existing fitted models/checkpoints.\n",
    "LUNG_RMST_BOOTSTRAPS (default 200), LUNG_RMST_WORKERS (default <=12), and LUNG_BENEFIT_OUTPUT are supported.\n", sep = "")
} else {
  stopifnot(file.exists("code/91_lung_pollution_survival_benefit.R"))
  Sys.setenv(LUNG_BENEFIT_HORIZON_DAYS = "365", LUNG_BENEFIT_BOOTSTRAPS = "0")
  stages <- c(if (!"--reuse-fitted" %in% args) c("code/00_prepare_lung_primary_inputs.R", "code/91_lung_pollution_survival_benefit.R"),
    "code/93_lung_pollution_rmst.R", "code/94_lung_pollution_rmst_bootstrap.R",
    "code/tests/validate_lung_survival_benefit_outputs.R", "code/tests/validate_lung_rmst_confidence_intervals.R")
  runner <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  for (stage in stages) {
    message(format(Sys.time(), "%F %T"), " | Primary analysis: ", stage)
    status <- system2(runner, stage)
    if (status != 0L) stop("Primary analysis stage failed: ", stage, " (exit ", status, ")")
  }
  message("Primary lung RMST workflow completed and validated.")
}
