source("code/revision_2026/common.R")
setDTthreads(2)
dest <- file.path(revision_dir, "las_era_lung")
cache <- file.path(revision_dir, "cache", "las_era_lung")
dir.create(dest, showWarnings = FALSE); dir.create(cache, showWarnings = FALSE)
cutoff <- as.Date("2023-03-08")
full <- baseline_data()
b <- copy(full[WL_ORG == "LU" & index_date <= cutoff])
original_n <- full[WL_ORG == "LU", .N]
b[, original_end_date := observed_end_date]
b[observed_end_date > cutoff, `:=`(adverse_event = 0L, transplant_or_improvement = 0L,
  other_exit = 0L, event_type = "administrative_censor")]
b[, `:=`(observed_end_date = pmin(observed_end_date, cutoff),
  multi_organ_candidate = as.integer(PERS_ID %in% full[, .(k = uniqueN(WL_ORG)), by = PERS_ID][k > 1, PERS_ID]))]
b[, followup_days := as.numeric(observed_end_date - index_date)]
rm(full); gc()
source("code/revision_2026/outcome_components.R")
components <- readRDS(file.path(revision_dir, "cache", "primary_outcome_components.rds"))
b <- merge(b, components[WL_ORG == "LU", .(PERS_ID, component)], by = "PERS_ID")
b[original_end_date > cutoff, component := "administrative_censor"]
b[, death_event := as.integer(component == "death")]
reg <- readRDS(file.path(revision_dir, "cache", "registrations.rds"))
reg <- merge(reg[WL_ORG == "LU"], b[, .(PERS_ID, cohort_start = index_date,
  cohort_end = observed_end_date, first_px = PX_ID)], by = "PERS_ID")
reg <- reg[index_date <= cutoff]
src <- readRDS(file.path(revision_dir, "cache", "official_lung_score_sources.rds"))
las <- unique(src$las[, .(PX_ID, score_date = CANHX_LAS_SCORE_BEGIN_DT,
  las = CANHX_LAS_SCORE_CALC)])
rm(src); gc()
las <- merge(las, reg[, .(PX_ID, PERS_ID, reg_start = index_date,
  reg_end = observed_end_date, cohort_start, cohort_end, first_px)], by = "PX_ID")
stopifnot(all(las$las[is.finite(las$las)] >= 0 & las$las[is.finite(las$las)] <= 100))
# Only scores from a registration already begun and not yet ended are eligible.
las <- las[score_date >= reg_start & score_date <= pmin(reg_end, cohort_end) & is.finite(las)]
las <- unique(las[, .(PERS_ID, PX_ID, first_px, score_date, las)])
# Without within-day timestamps, conflicting values are not arbitrarily ordered.
conflicts <- las[, .N, by = .(PX_ID, score_date)][N > 1]
las <- las[!conflicts, on = .(PX_ID, score_date)]
# For simultaneous center updates, prefer the cohort's first registration on ties.
las[, index_registration := PX_ID == first_px]
setorder(las, PERS_ID, score_date, -index_registration, PX_ID)
las <- unique(las, by = c("PERS_ID", "score_date"))
setorder(las, PERS_ID, score_date)
las[, changed := is.na(shift(las)) | las != shift(las), by = PERS_ID]
events <- las[changed == TRUE, .(PERS_ID, score_date, las)]
setkey(events, PERS_ID, score_date)
b[, baseline_las := events[.SD, on = .(PERS_ID, score_date = index_date), roll = Inf, x.las]]
b[, baseline_las_date := events[.SD, on = .(PERS_ID, score_date = index_date), roll = Inf, x.score_date]]
b[age < 12, baseline_las := NA_real_]
ref <- median(b[is.finite(baseline_las) & baseline_las > 0, baseline_las])
stopifnot(is.finite(ref), all(b$baseline_las_date <= b$index_date, na.rm = TRUE))
fwrite(data.table(measure = c("Original lung candidates", "Excluded listings after March 8 2023",
  "Retained lung candidates", "Follow-up truncated at transition", "Adverse outcomes", "Death outcomes",
  "Conflicting within-registration score dates omitted", "Eligible score-change records",
  "Baseline LAS unavailable or age under 12", "Baseline recorded zero LAS", "Baseline age under 12",
  "LAS imputation reference"), value = c(original_n, original_n - nrow(b), nrow(b),
  sum(b$original_end_date > cutoff), sum(b$adverse_event), sum(b$death_event), nrow(conflicts), nrow(events),
  sum(!is.finite(b$baseline_las)), sum(b$baseline_las == 0, na.rm = TRUE), sum(b$age < 12), ref)),
  file.path(dest, "cohort_and_las_audit.csv"))
summary <- b[, .(N = .N, observed = sum(is.finite(baseline_las)),
  median = median(baseline_las, na.rm = TRUE), q25 = quantile(baseline_las, .25, na.rm = TRUE),
  q75 = quantile(baseline_las, .75, na.rm = TRUE))]
print(summary); fwrite(summary, file.path(dest, "baseline_las_summary.csv"))
fwrite(b[, .N, by = component], file.path(dest, "outcome_components.csv"))
saveRDS(list(cohort = b, scores = events, imputation = ref, cutoff = cutoff),
  file.path(cache, "prepared.rds"), compress = FALSE)
print(fread(file.path(dest, "cohort_and_las_audit.csv")))
