# Publication cohort: retain original nonlung rows and the audited LAS-era lung cohort.
read_publication_cohort <- function(path) {
  dat <- readr::read_csv(path, show_col_types = FALSE)
  prep <- readRDS("output/revision_2026_checklist/cache/las_era_lung/prepared.rds")
  lung <- tibble::as_tibble(prep$cohort)
  stopifnot(all(lung$index_date <= prep$cutoff), all(lung$observed_end_date <= prep$cutoff))
  # Retain source column types, especially categorical fields and residential ZIP.
  dat <- dat[dat$WL_ORG != "LU" | dat$PERS_ID %in% lung$PERS_ID, ]
  i <- which(dat$WL_ORG == "LU")
  j <- match(dat$PERS_ID[i], lung$PERS_ID)
  for (nm in c("observed_end_date", "adverse_event", "transplant_or_improvement", "other_exit", "event_type", "followup_days")) {
    if (nm %in% names(dat)) dat[[nm]][i] <- lung[[nm]][j]
  }
  dat$baseline_las <- NA_real_
  dat$baseline_las[i] <- lung$baseline_las[j]
  dat$lung_las_era <- dat$WL_ORG == "LU"
  dat$candidate_zip <- sprintf("%05d", as.integer(dat$candidate_zip))
  stopifnot(length(i) == 42874L, nrow(dat) == 1007801L)
  dat
}

publication_baseline_data <- function() {
  old <- readRDS("output/revision_2026_checklist/cache/baseline.rds")
  lung <- readRDS("output/revision_2026_checklist/cache/las_era_lung/prepared.rds")$cohort
  data.table::rbindlist(list(old[WL_ORG != "LU"], lung), fill = TRUE)
}

publication_period_note <- paste(
  "Lung listings and follow-up end March 8, 2023; other organs retain their original study periods.",
  "Lung time-varying models use official calculated LAS updated at monthly interval starts; CAS is not used.")
