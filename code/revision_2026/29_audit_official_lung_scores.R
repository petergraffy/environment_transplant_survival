source("code/revision_2026/common.R")
setDTthreads(2)
dest <- file.path(revision_dir, "official_lung_scores")
dir.create(dest, showWarnings = FALSE)
path <- "C:/Users/Peter Graffy/Box/Q1 2026 SAF/Lung_supplement2603"
b <- baseline_data()[WL_ORG == "LU", .(PERS_ID, PX_ID, index_date, observed_end_date, episode_px_ids, event_type)]
cas <- as.data.table(haven::read_sas(file.path(path, "cas_2603.sas7bdat")))
setnames(cas, c("pers_id", "px_id"), c("PERS_ID", "PX_ID"))
summary <- cas[, .(rows = .N, people = uniqueN(PERS_ID), registrations = uniqueN(PX_ID),
  dates = uniqueN(CAS_dt), min_date = min(CAS_dt, na.rm = TRUE), max_date = max(CAS_dt, na.rm = TRUE),
  missing_score = sum(is.na(CAS)), min_score = min(CAS, na.rm = TRUE),
  median_score = median(CAS, na.rm = TRUE), max_score = max(CAS, na.rm = TRUE),
  distance_values = uniqueN(DistanceInNM), distance_missing = sum(is.na(DistanceInNM)),
  distance_min = min(DistanceInNM, na.rm = TRUE), distance_max = max(DistanceInNM, na.rm = TRUE))]
print(summary); fwrite(summary, file.path(dest, "cas_source_summary.csv"))
fwrite(cas[, .N, by = CAS_dt], file.path(dest, "cas_dates.csv"))
j <- merge(cas, b[, .(PERS_ID, index_date, observed_end_date)], by = "PERS_ID")
coverage <- j[, .(linked_people = uniqueN(PERS_ID), linked_rows = .N,
  people_at_or_before_listing = uniqueN(PERS_ID[CAS_dt <= index_date]),
  people_during_followup = uniqueN(PERS_ID[CAS_dt >= index_date & CAS_dt <= observed_end_date]),
  people_after_followup = uniqueN(PERS_ID[CAS_dt > observed_end_date]))]
print(coverage); fwrite(coverage, file.path(dest, "cas_cohort_timing.csv"))
log_revision("Reading LAS score history")
las <- as.data.table(haven::read_sas(file.path(path, "las_saf2603.sas7bdat")))
s <- las[, .(rows = .N, registrations = uniqueN(PX_ID),
  min_date = min(CANHX_LAS_SCORE_BEGIN_DT, na.rm = TRUE), max_date = max(CANHX_LAS_SCORE_BEGIN_DT, na.rm = TRUE),
  missing_date = sum(is.na(CANHX_LAS_SCORE_BEGIN_DT)),
  missing_calculated = sum(is.na(CANHX_LAS_SCORE_CALC)), missing_match = sum(is.na(CANHX_LAS_SCORE_MATCH)),
  min_calculated = min(CANHX_LAS_SCORE_CALC, na.rm = TRUE), max_calculated = max(CANHX_LAS_SCORE_CALC, na.rm = TRUE))]
print(s); fwrite(s, file.path(dest, "las_source_summary.csv"))
saveRDS(list(las = las, cas = cas), file.path(revision_dir, "cache", "official_lung_score_sources.rds"), compress = FALSE)
yearly <- las[, .(rows = .N, registrations = uniqueN(PX_ID),
  observed_calculated = sum(!is.na(CANHX_LAS_SCORE_CALC)),
  median_calculated = median(CANHX_LAS_SCORE_CALC, na.rm = TRUE)),
  by = .(year = as.integer(format(CANHX_LAS_SCORE_BEGIN_DT, "%Y")))]
setorder(yearly, year); print(yearly)
fwrite(yearly, file.path(dest, "las_coverage_by_year.csv"))
source("code/saf_paths.R")
sp <- get_saf_paths(release = "q1_2026")
regs <- as.data.table(haven::read_sas(file.path(sp$pubsaf_dir, "cand_thor.sas7bdat"),
  col_select = tidyselect::all_of(c("PERS_ID", "PX_ID", "WL_ORG", "CAN_LISTING_DT", "CAN_REM_DT", "CAN_REM_CD"))))
regs <- unique(regs[WL_ORG == "LU"])
stopifnot(!anyDuplicated(regs[, .(PERS_ID, PX_ID)]))
timing <- merge(cas, regs, by = c("PERS_ID", "PX_ID"))
timing[, `:=`(days_after_listing = as.numeric(CAS_dt - CAN_LISTING_DT), days_before_removal = as.numeric(CAN_REM_DT - CAS_dt))]
timing_summary <- timing[, .(N = .N, at_listing = sum(days_after_listing == 0, na.rm = TRUE),
  before_listing = sum(days_after_listing < 0, na.rm = TRUE), after_listing = sum(days_after_listing > 0, na.rm = TRUE),
  same_as_removal = sum(days_before_removal == 0, na.rm = TRUE),
  after_removal = sum(days_before_removal < 0, na.rm = TRUE),
  days_after_listing_median = median(days_after_listing, na.rm = TRUE),
  days_before_removal_median = median(days_before_removal, na.rm = TRUE)), by = CAN_REM_CD]
print(timing_summary); fwrite(timing_summary, file.path(dest, "cas_registration_date_audit.csv"))
available <- unique(j[CAS_dt >= index_date & CAS_dt <= observed_end_date, PERS_ID])
b[, cas_during_followup := PERS_ID %in% available]
outcome_coverage <- b[, .(N = .N, CAS_available = sum(cas_during_followup),
  percent = 100 * mean(cas_during_followup)), by = event_type]
fwrite(outcome_coverage, file.path(dest, "cas_coverage_by_primary_outcome.csv"))
