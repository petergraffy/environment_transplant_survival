source("code/revision_2026/common.R")
source("code/publication_cohort.R")
setDTthreads(2)
dest <- file.path(revision_dir, "las_era_lung")
backup <- file.path(revision_dir, "cache", "las_publication_before")
dir.create(backup, recursive = TRUE, showWarnings = FALSE)
original <- function(path) {
  p <- file.path(backup, gsub("[/\\\\]", "__", path))
  if (!file.exists(p)) file.copy(path, p)
  fread(p)
}
replace_lung <- function(path, rows, key = "organ", lung = "LU") {
  old <- original(path)
  rows <- as.data.table(rows)[, intersect(names(rows), union(names(old), "warnings")), with = FALSE]
  out <- rbindlist(list(old[get(key) != lung], rows), use.names = TRUE, fill = TRUE)
  before <- old[get(key) != lung]
  after <- out[get(key) != lung]
  stopifnot(all(vapply(names(old), function(nm)
    identical(as.character(before[[nm]]), as.character(after[[nm]])), logical(1))))
  fwrite(out, path)
}
r <- fread(file.path(dest, "pollutant_results.csv"))
r[, pollutant := sub("_.*", "", term)]
convert <- function(z) {
  z[, .(organ = "LU", organ_label = "Lung", exposure = pollutant, pollutant,
    n = N, people = N, candidate_organ_episodes = N, intervals,
    adverse_events = events, hazard_ratio = HR, conf_low = HR_low, conf_high = HR_high,
    p_value = p.value, adjustment_set = formula, warnings = fifelse(status == "OK", "", status))]
}
for (m in c("Baseline", "Time-varying")) {
  path <- if(m == "Baseline") "output/prior_year_pollution_cox_svi/prior_year_pollution_cox_svi_results.csv" else
    "output/timevarying_pollution_cox_svi/timevarying_pollution_cox_svi_results.csv"
  z <- convert(r[model == m & variant == "Primary"])
  z[, `:=`(variance_estimator = "model_based", center_adjustment = "stratified_baseline_hazard_by_listing_center")]
  replace_lung(path, z)
}
sub <- r[variant %in% c("Age subgroup", "Sex subgroup", "Race subgroup")]
z <- convert(sub)
z[, `:=`(model = fifelse(sub$model == "Baseline", "baseline_prior_year", "time_varying"),
  subgroup_label = c("Age subgroup"="Age group", "Sex subgroup"="Sex", "Race subgroup"="Race")[sub$variant],
  subgroup_level = sub$subgroup, error = NA_character_,
  hazard_ratio_ci = sprintf("%.2f (%.2f-%.2f)", hazard_ratio, conf_low, conf_high))]
replace_lung("output/prior_year_and_timevarying_pollution_subgroup_cox/prior_year_and_timevarying_pollution_subgroup_cox_results.csv", z)
variants <- c("No center strata"="Single-pollutant model without listing center",
  "Single organ only"="Single-pollutant model restricted to single-organ candidates",
  "Adjust multi-organ"="Single-pollutant model additionally adjusted for multi-organ candidate status",
  "PM2.5 + NO2"="Multipollutant model: PM2.5 + NO2", "PM2.5 + NO2 + O3"="Multipollutant model: PM2.5 + NO2 + O3")
s <- r[variant %in% names(variants)]
z <- convert(s)
z[, `:=`(analysis = paste(s$model, "cause-specific Cox"), model = unname(variants[s$variant]),
  exposure_label = c(pm25="PM2.5", no2="NO2", o3="O3")[s$pollutant])]
replace_lung("output/updated_methods_sensitivity_models/updated_methods_sensitivity_results.csv", z)

# Keep matched death-only comparisons synchronized with the new primary models.
for (filename in c("pollutant_results.csv", "full_model_results.csv")) {
  x <- if(filename == "pollutant_results.csv") copy(r) else fread(file.path(dest, "full_coefficients.csv"))
  x <- x[variant %in% c("Primary", "Death only")]
  x[, `:=`(organ = "LU", endpoint = fifelse(variant == "Primary", "Composite", "Death only"),
    candidate_organ_episodes = N, events = events, hazard_ratio = HR, warnings = fifelse(status == "OK", "", status))]
  if(filename == "full_model_results.csv") {
    # Each single-pollutant formula identifies its pollutant, including nuisance terms.
    x[, pollutant := fifelse(grepl("pm25_", formula), "pm25", fifelse(grepl("no2_", formula), "no2", "o3"))]
  }
  replace_lung(file.path(revision_dir, "death_only_sensitivity", filename), x)
  if(filename == "full_model_results.csv") for(m in unique(x$model)) for(p in unique(x$pollutant)) {
    fwrite(x[model==m & pollutant==p],file.path(revision_dir,"death_only_sensitivity",paste0(tolower(m),"_LU_",p,"_coefficients.csv")))
  }
}
dat <- publication_baseline_data()
fwrite(dat[, .(people=.N), by=WL_ORG], file.path(dest, "publication_cohort_counts.csv"))

# Freeze national quartiles to the reference cohort; update lung curves only.
ref <- baseline_data()
raws <- c("PM2.5"="pm25_prior_ug_m3", "NO2"="no2_prior_ppb", "O3"="o3_prior_ppb")
quartiles <- c("Q1 lowest", "Q2", "Q3", "Q4 highest")
curves <- list(); times <- list(); summaries <- list(); cutpoints <- list()
for(p in names(raws)) {
  col <- raws[[p]]
  eligible_ref <- ref[is.finite(get(col)) & followup_days > 0 & !is.na(adverse_event) & !is.na(transplant_or_improvement)]
  cuts <- as.numeric(quantile(eligible_ref[[col]], seq(0,1,.25))); cuts[c(1,5)] <- c(-Inf,Inf)
  cutpoints[[p]] <- data.table(pollutant=p, q25=cuts[2], q50=cuts[3], q75=cuts[4])
  d <- copy(dat[WL_ORG == "LU" & is.finite(get(col)) & followup_days > 0])
  d[, quartile := cut(get(col), cuts, labels=quartiles)]
  d[, event := factor(fcase(adverse_event==1,"adverse",transplant_or_improvement==1,"transplant_or_improvement",default="censor"),
    levels=c("censor","adverse","transplant_or_improvement"))]
  fit <- survfit(Surv(followup_days/365.25, event) ~ quartile, data=as.data.frame(d))
  tidy <- function(s) {
    data.table(time=s$time, cif_adverse=s$pstate[,match("adverse",s$states)],
      cif_transplant_or_improvement=s$pstate[,match("transplant_or_improvement",s$states)],
      n_risk=if(is.matrix(s$n.risk)) s$n.risk[,1] else s$n.risk,
      strata=as.character(s$strata), quartile=sub("^quartile=","",as.character(s$strata)),
      pollutant=p, unit_label=if(p=="PM2.5") "ug/m3" else "ppb", organ="Lung")
  }
  curves[[p]] <- tidy(summary(fit))
  times[[p]] <- tidy(summary(fit,times=c(1,3,5,10),extend=TRUE))
  times[[p]][, `:=`(cif_adverse_percent=100*cif_adverse, cif_transplant_or_improvement_percent=100*cif_transplant_or_improvement)]
  summaries[[p]] <- d[, .(n=.N,adverse_events=sum(adverse_event),competing_events=sum(transplant_or_improvement),
    censored=sum(event=="censor"),median_prior_exposure=median(get(col))),by=quartile][,
      `:=`(pollutant=p,unit_label=if(p=="PM2.5") "ug/m3" else "ppb",organ="Lung")]
}
aj <- "output/prior_year_pollution_quartile_aalen_johansen_cif"
replace_lung(file.path(aj,"prior_year_pollution_quartile_aalen_johansen_curve_data.csv"),rbindlist(curves),"organ","Lung")
replace_lung(file.path(aj,"prior_year_pollution_quartile_cif_at_1_3_5_10_years.csv"),rbindlist(times),"organ","Lung")
replace_lung(file.path(aj,"prior_year_pollution_quartile_event_summary.csv"),rbindlist(summaries),"organ","Lung")
fwrite(rbindlist(cutpoints),file.path(dest,"national_reference_quartile_cutpoints.csv"))

# Map weights refer to candidates' first listing center, as in the primary models.
path <- "output/center_maps/transplant_centers_by_organ_mapped.csv"
old <- original(path)
counts <- dat[WL_ORG=="LU",.(n_candidates=.N),by=.(center_code=as.character(listing_center))]
z <- merge(old[organ=="LU", !"n_candidates"], counts, by="center_code", all=FALSE)
replace_lung(path,z)
all_centers <- fread(path)
path2 <- "output/center_maps/transplant_centers_mapped.csv"
old2 <- original(path2)
totals <- all_centers[, .(total_candidates=sum(n_candidates),dominant_organ=organ_label[which.max(n_candidates)],n_organs=.N),by=center_code]
fwrite(merge(old2[,!c("total_candidates","dominant_organ","n_organs")],totals,by="center_code"),path2)
writeLines(c(publication_period_note,
  "National AJ quartile thresholds are frozen to the original full reference cohort.",
  "Nonlung source result rows were checked for exact equality before export."), file.path(dest,"publication_input_validation.txt"))
log_revision("Publication model inputs and AJ lung estimates synchronized")
