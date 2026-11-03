source("code/revision_2026/common.R")
source("code/saf_paths.R")
setDTthreads(2)
dest <- file.path(revision_dir, "retransplant_candidates")
dir.create(dest, showWarnings = FALSE)
b <- baseline_data()[, .(PERS_ID, PX_ID, WL_ORG, index_date)]
paths <- get_saf_paths(release = "q1_2026")
fields <- c("CAN_PREV_HR", "CAN_PREV_HL", "CAN_PREV_KI", "CAN_PREV_KP",
            "CAN_PREV_LI", "CAN_PREV_LU", "CAN_PREV_PA", "CAN_PREV_IN", "CAN_PREV_TX")
raw <- rbindlist(lapply(c("cand_kipa", "cand_liin", "cand_thor"), function(nm) {
  as.data.table(haven::read_sas(file.path(paths$pubsaf_dir, paste0(nm, ".sas7bdat")),
    col_select = tidyselect::any_of(c("PERS_ID", "PX_ID", "WL_ORG", fields))))
}), fill = TRUE)
raw <- unique(raw[WL_ORG %in% names(organ_names)])
stopifnot(!anyDuplicated(raw[, .(PERS_ID, PX_ID, WL_ORG)]))
d <- merge(b, raw, by = c("PERS_ID", "PX_ID", "WL_ORG"), all.x = TRUE)
stopifnot(nrow(d) == nrow(b))
freq <- rbindlist(lapply(fields, function(v) d[, .N, by = c("WL_ORG", v)][, .(organ = WL_ORG, field = v, value = as.numeric(get(v)), N)]))
fwrite(freq, file.path(dest, "prior_transplant_field_frequencies.csv"))
stopifnot(all(freq$value[!is.na(freq$value)] %in% c(0, 1)))
print(freq[!is.na(value) & value != 0])
same <- list(HR = c("CAN_PREV_HR", "CAN_PREV_HL"), KI = c("CAN_PREV_KI", "CAN_PREV_KP"),
             LI = "CAN_PREV_LI", LU = c("CAN_PREV_LU", "CAN_PREV_HL"))
for (org in names(same)) {
  cols <- same[[org]]
  d[WL_ORG == org, `:=`(same_organ_prior = rowSums(.SD > 0, na.rm = TRUE) > 0,
    same_organ_fields_missing = rowSums(!is.na(.SD)) == 0), .SDcols = cols]
}
organ_fields <- setdiff(fields, "CAN_PREV_TX")
d[, any_organ_prior := rowSums(.SD > 0, na.rm = TRUE) > 0, .SDcols = organ_fields]
out <- d[, .(N = .N, same_organ_prior_n = sum(same_organ_prior),
  same_organ_prior_percent = 100 * mean(same_organ_prior),
  same_organ_fields_all_missing = sum(same_organ_fields_missing),
  any_organ_prior_n = sum(any_organ_prior), previous_transplants_field_positive = sum(CAN_PREV_TX > 0, na.rm = TRUE)), by = WL_ORG]
fwrite(out, file.path(dest, "retransplant_candidates_by_organ.csv"))
print(out)
overall <- d[, .(candidate_organ_observations = .N, same_organ_prior_n = sum(same_organ_prior),
  same_organ_prior_percent = 100 * mean(same_organ_prior),
  unique_people_with_same_organ_prior = uniqueN(PERS_ID[same_organ_prior]),
  any_organ_prior_n = sum(any_organ_prior))]
fwrite(overall, file.path(dest, "retransplant_candidates_overall.csv"))
print(overall)
writeLines(c("# Retransplant candidates at baseline", "",
  "Counts use the first listing PX_ID retained in the primary deduplicated candidate-organ cohort, not later registration records or subsequent transplant outcomes.",
  "Same-organ prior transplant: heart CAN_PREV_HR or CAN_PREV_HL; kidney CAN_PREV_KI or CAN_PREV_KP; liver CAN_PREV_LI; lung CAN_PREV_LU or CAN_PREV_HL. Fields are binary 0/1 in this cohort.",
  "Any-organ history uses the union of CAN_PREV_HR/HL/KI/KP/LI/LU/PA/IN. Counts agree with CAN_PREV_TX.",
  "These are documented baseline histories in the candidate files, not counts of repeat transplant procedures during follow-up. They have not been independently validated against longitudinal recipient/transplant histories.",
  "Aggregate totals count candidate-organ observations; unique people are also reported separately. No primary analyses or Table 1 were changed."),
  file.path(dest, "README.md"))
