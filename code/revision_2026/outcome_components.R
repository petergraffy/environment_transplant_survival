# Resolve the same terminal registration as the primary candidate-level cohort.
primary_outcome_components <- function(base, cache_dir = "output/revision_2026_checklist/cache") {
  suppressPackageStartupMessages(library(data.table))
  b <- as.data.table(base)
  cache_file <- file.path(cache_dir, "primary_outcome_components.rds")
  if (file.exists(cache_file)) {
    out <- readRDS(cache_file)
  } else {
    registration_file <- file.path(cache_dir, "registrations.rds")
    if (!file.exists(registration_file)) stop(
      "Run code/revision_2026/02_prepare_registration_sensitivities.R first to prepare the private registration cache.")
    source("code/saf_paths.R")
    paths <- get_saf_paths(release = "q1_2026")
    reg <- copy(readRDS(registration_file))
    raw <- rbindlist(lapply(c("cand_kipa", "cand_liin", "cand_thor"), function(nm) {
      as.data.table(haven::read_sas(file.path(paths$pubsaf_dir, paste0(nm, ".sas7bdat")),
        col_select = tidyselect::all_of(c("PERS_ID", "PX_ID", "WL_ORG", "CAN_REM_CD"))))
    }))
    raw <- unique(raw[WL_ORG %in% c("HR", "KI", "LI", "LU")])
    stopifnot(!anyDuplicated(raw[, .(PERS_ID, PX_ID, WL_ORG)]))
    reg <- merge(reg, raw, by = c("PERS_ID", "PX_ID", "WL_ORG"), all.x = TRUE)
    reg[, priority := fifelse(adverse_event == 1, 3L,
      fifelse(transplant_or_improvement == 1, 2L, fifelse(other_exit == 1, 1L, 0L)))]
    reg[, has_terminal := any(priority > 0), by = .(PERS_ID, WL_ORG)]
    end <- reg[priority > 0 | !has_terminal]
    setorder(end, PERS_ID, WL_ORG, observed_end_date, priority, waitlist_row_id)
    end <- end[end[, .I[.N], by = .(PERS_ID, WL_ORG)]$V1]
    end[, component := fcase(
      priority == 3 & CAN_REM_CD == 8, "death",
      priority == 3 & CAN_REM_CD == 13, "deterioration",
      priority == 2 & CAN_REM_CD == 12, "improvement",
      priority == 2, "transplant",
      priority == 1, "other_exit", default = "administrative_censor")]
    out <- end[, .(PERS_ID, WL_ORG, observed_end_date, adverse_event,
      transplant_or_improvement, other_exit, terminal_removal_code = CAN_REM_CD,
      component, death_event = as.integer(component == "death"))]
  }
  # Cache the original cohort resolution, never the publication-specific censoring.
  if (!file.exists(cache_file)) saveRDS(out, cache_file, compress = FALSE)
  out <- copy(out)
  if ("lung_las_era" %in% names(b) && any(b$lung_las_era)) {
    cutoff <- as.Date("2023-03-08")
    out <- out[WL_ORG != "LU" | PERS_ID %in% b[WL_ORG == "LU", PERS_ID]]
    out[WL_ORG == "LU" & observed_end_date > cutoff, `:=`(
      observed_end_date = cutoff, adverse_event = 0L, transplant_or_improvement = 0L,
      other_exit = 0L, death_event = 0L, component = "administrative_censor")]
  }
  check <- merge(b[, .(PERS_ID, WL_ORG, observed_end_date, adverse_event,
    transplant_or_improvement, other_exit)], out, by = c("PERS_ID", "WL_ORG"),
    suffixes = c("_primary", "_resolved"))
  stopifnot(nrow(check) == nrow(b), !anyDuplicated(check[, .(PERS_ID, WL_ORG)]))
  for (v in c("observed_end_date", "adverse_event", "transplant_or_improvement", "other_exit")) {
    stopifnot(all(check[[paste0(v, "_primary")]] == check[[paste0(v, "_resolved")]]))
  }
  stopifnot(all(out$adverse_event == as.integer(out$component %in% c("death", "deterioration"))),
    all(out$transplant_or_improvement == as.integer(out$component %in% c("transplant", "improvement"))))
  out[]
}
