#!/usr/bin/env Rscript
# Adapter for the established cohort definitions, not a new cohort specification.
# Build only the private inputs needed by the primary lung workflow; do not fit
# the historical four-organ adverse-outcome models or registration sensitivities.
source("code/r_runtime.R"); ensure_user_library()
suppressPackageStartupMessages({library(data.table); library(haven)})
source("code/saf_paths.R")
archive <- "code/archive/pre_lung_rmst"
source(file.path(archive, "revision_2026", "common.R"))
setDTthreads(2)
prepared_path <- "output/revision_2026_checklist/cache/las_era_lung/prepared.rds"
if (file.exists(prepared_path)) {
  prepared <- readRDS(prepared_path)
  stopifnot(all(c("cohort", "scores", "imputation", "cutoff") %in% names(prepared)),
    !anyDuplicated(prepared$cohort$PERS_ID), all(prepared$cohort$WL_ORG == "LU"),
    all(prepared$cohort$index_date <= as.Date("2023-03-08")),
    all(prepared$cohort$observed_end_date <= as.Date("2023-03-08")),
    is.finite(prepared$imputation), as.Date(prepared$cutoff) == as.Date("2023-03-08"))
  message("Verified established LAS-era cohort cache: ", nrow(prepared$cohort), " candidates.")
} else {
  analysis_path <- "output/primary_waitlist_period_pollution_cox/primary_waitlist_period_pollution_analysis_dataset.csv.gz"
  if (!file.exists(analysis_path)) {
    log_revision("Preparing established deduplicated cohort without historical model fits")
    env <- load_script_prefix(file.path(archive, "50_primary_waitlist_period_pollution_cox.R"), "result_paths")
    rm(env); gc(FALSE)
  }
  base <- baseline_data()
  registration_path <- file.path(revision_dir, "cache", "registrations.rds")
  if (!file.exists(registration_path)) {
    env <- load_script_prefix(file.path(archive, "50_primary_waitlist_period_pollution_cox.R"), "registration_cohort")
    cols <- c("PERS_ID", "PX_ID", "WL_ORG", "waitlist_row_id", "index_date", "observed_end_date",
      "candidate_zip", "age", "sex", "race", "listing_center", "listing_year", "adverse_event",
      "transplant_or_improvement", "other_exit")
    reg <- as.data.table(env$cohort)[, ..cols]
    stopifnot(!anyDuplicated(reg[, .(PERS_ID, PX_ID, WL_ORG)]))
    setorder(reg, PERS_ID, WL_ORG, index_date, observed_end_date, waitlist_row_id)
    reg[, spell := spell_ids(index_date, observed_end_date), by = .(PERS_ID, WL_ORG)]
    saveRDS(reg, registration_path, compress = FALSE)
    rm(env, reg); gc(FALSE)
  }
  source(file.path(archive, "revision_2026", "outcome_components.R"))
  invisible(primary_outcome_components(base))
  score_path <- file.path(revision_dir, "cache", "official_lung_score_sources.rds")
  if (!file.exists(score_path)) {
    paths <- get_saf_paths(release = "q1_2026")
    supplement <- file.path(paths$saf_dir, "Lung_supplement2603")
    las <- as.data.table(read_sas(file.path(supplement, "las_saf2603.sas7bdat")))
    # Preserve the shared source-cache shape for archived audits. CAS is not
    # modeled in the current analysis.
    cas_path <- file.path(supplement, "cas_2603.sas7bdat")
    cas <- if (file.exists(cas_path)) as.data.table(read_sas(cas_path)) else NULL
    if (!is.null(cas) && all(c("pers_id", "px_id") %in% names(cas)))
      setnames(cas, c("pers_id", "px_id"), c("PERS_ID", "PX_ID"))
    saveRDS(list(las = las, cas = cas), score_path, compress = FALSE)
    rm(las, cas); gc(FALSE)
  }
  options(las.bundle = "las_era_lung", las.extend_followup = FALSE)
  source(file.path(archive, "revision_2026", "30_prepare_las_era_lung.R"))
  stopifnot(file.exists(prepared_path))
}
