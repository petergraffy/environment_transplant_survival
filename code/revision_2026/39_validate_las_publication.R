source("code/revision_2026/common.R")
source("code/publication_cohort.R")
library(gt)
dest <- file.path(revision_dir,"las_era_lung")
cohort <- fread(file.path(dest,"publication_cohort_counts.csv"))
stopifnot(sum(cohort$people)==1007801L,cohort[WL_ORG=="LU",people]==42874L)
outcomes <- fread(file.path(dest,"publication_outcomes_by_organ.csv"))
components <- c("death","deterioration","transplant","improvement","other_exit","administrative_censor")
stopifnot(all(rowSums(outcomes[,..components])==outcomes$N),outcomes[WL_ORG=="LU",death+deterioration]==5545L)
base <- fread("output/prior_year_pollution_cox_svi/prior_year_pollution_cox_svi_results.csv")
tv <- fread("output/timevarying_pollution_cox_svi/timevarying_pollution_cox_svi_results.csv")
new <- fread(file.path(dest,"pollutant_results.csv"))[variant=="Primary"]
new[,pollutant:=sub("_.*","",term)]
for(p in c("pm25","no2","o3")) {
  stopifnot(abs(base[organ=="LU" & exposure==p,hazard_ratio]-new[model=="Baseline" & pollutant==p,HR])<1e-12,
    abs(tv[organ=="LU" & pollutant==p,hazard_ratio]-new[model=="Time-varying" & pollutant==p,HR])<1e-12)
}
bd <- fread(file.path(revision_dir,"baseline_diagnostics","baseline_reproduction.csv"))
td <- fread(file.path(revision_dir,"timevarying_diagnostics","timevarying_reproduction.csv"))
check <- merge(base,bd,by.x=c("organ","exposure"),by.y=c("organ","pollutant"))
stopifnot(nrow(check)==12,all(abs(check$hazard_ratio.x-check$hazard_ratio.y)<1e-7),all(check$n.x==check$n.y))
check <- merge(tv,td,by=c("organ","pollutant"))
stopifnot(nrow(check)==12,all(abs(check$hazard_ratio-check$estimate)<1e-7),all(check$candidate_organ_episodes==check$n))
aj <- fread(file.path(revision_dir,"aj_with_risk_tables","cohort_time_zero_reconciliation.csv"))
stopifnot(all(aj$aj_n==aj$time_zero_total),all(aj$cohort_n==aj$aj_n+aj$excluded_followup+aj$excluded_outcome+aj$excluded_exposure))
path <- "output/tables/table1_baseline_waitlist_characteristics/table1_baseline_characteristics_by_organ.csv"
current <- fread(path)
original_text <- system2("git",c("show",paste0("HEAD:",path)),stdout=TRUE)
old <- fread(text=paste(original_text,collapse="\n"))
matched <- merge(current,old,by=c("section","characteristic"),suffixes=c("_new","_old"))
for(org in c("Heart","Kidney","Liver")) stopifnot(identical(matched[[paste0(org,"_new")]],matched[[paste0(org,"_old")]]))
stopifnot(current[grepl("Lung allocation score",characteristic),Lung]=="37.7 (33.9 to 45.6)")
maps <- fread("output/figures/pollution_absolute_1yr_survival_maps/pollution_absolute_1yr_survival_map_summary.csv")
stopifnot(all(maps[organ=="Lung",model_n]==40721),all(maps[organ=="Lung",adverse_events]==5206))
candidate_map <- fread("output/figures/candidates_heart_lung_absolute_risk_abc/candidates_heart_lung_absolute_risk_abc_summary.csv")
stopifnot(candidate_map$total_waitlist_episodes==1007801L)

# Copyable HTML versions of the numeric sensitivity tables.
for(path in c("output/tables/updated_methods_sensitivity_models/supplemental_baseline_sensitivity_wide_table.csv",
  "output/tables/updated_methods_sensitivity_models/supplemental_timevarying_sensitivity_wide_table.csv",
  "output/prior_year_pollution_baseline_severity_associations/prior_year_pollution_baseline_severity_associations_table.csv",
  "output/timevarying_pollution_severity_models_with_listing_year/timevarying_pollution_severity_table.csv")) {
  gtsave(gt(as.data.frame(fread(path))) |> tab_source_note(publication_period_note) |>
    tab_options(table.font.names="Arial",table.font.size=px(13)),sub("[.]csv$",".html",path))
}

entries <- list()
add <- function(group,folder,pattern,status="Updated") {
  files <- list.files(folder,pattern=pattern,full.names=TRUE)
  if(group=="Candidate and thoracic risk ABC") files <- files[!grepl("_updated[.]pdf$",files)]
  stopifnot(length(files)>0,all(file.info(files)$size>0))
  entries[[length(entries)+1L]] <<- data.table(group=group,path=files,status=status,bytes=file.info(files)$size)
}
add("Main Cox forest","output/figures/figure1_cox_baseline_timevarying","[.](pdf|png)$")
add("AJ curves with risk tables",file.path(revision_dir,"aj_with_risk_tables"),"[.](pdf|png|csv)$")
add("Simplified AJ curves","output/figures/prior_year_pollution_quartile_aalen_johansen_simplified","[.](pdf|png)$")
add("Subgroup forests","output/figures/prior_year_and_timevarying_pollution_subgroup_cox/las_era_20261015","[.](pdf|png)$")
add("Cohort flow","output/figures/cohort_flow_diagram","[.](pdf|png|csv|txt)$")
add("Absolute-scale maps","output/figures/pollution_absolute_1yr_survival_maps","[.](pdf|png|csv)$")
add("Candidate and thoracic risk ABC","output/figures/candidates_heart_lung_absolute_risk_abc","[.](pdf|png|csv)$")
add("Center maps","output/figures/transplant_centers_by_organ_four_panel","[.](pdf|png|csv)$")
add("Candidate maps","output/figures/waitlist_candidates_zcta_map","[.](pdf|png|csv)$")
add("Candidate-center AB","output/figures/waitlist_geography_ab_maps","[.](pdf|png|csv)$")
add("Table 1","output/tables/table1_baseline_waitlist_characteristics","[.](html|csv|txt)$")
add("Covariate missingness","output/tables/primary_cox_model_variables_missingness","[.](html|csv|txt)$")
add("Score definitions","output/tables/organ_score_formula_weights","[.](html|csv)$")
add("Sensitivity models","output/tables/updated_methods_sensitivity_models","[.](html|csv|tsv|txt)$")
add("Three-tier baseline models","output/tables/prior_year_initial_listing_cox_model_set","[.](csv|tsv|txt)$")
add("Three-tier baseline forests","output/figures/prior_year_initial_listing_cox_model_set","[.](pdf|png)$")
add("Death-only sensitivity",file.path(revision_dir,"death_only_sensitivity"),"(comparison[.]csv|results[.]csv|[.]html)$")
add("Reviewer supplement",file.path(revision_dir,"supplement"),"[.](html|csv)$")
add("PH diagnostic figures",file.path(revision_dir,"baseline_diagnostics"),"[.](pdf|png)$")
add("PH diagnostic figures",file.path(revision_dir,"timevarying_diagnostics"),"[.](pdf|png)$")
add("Severity-response figure","output/figures/pollution_severity_associations","^figure4_.*[.](pdf|png)$")
add("Severity-response tables","output/prior_year_pollution_baseline_severity_associations","[.](html|csv)$")
add("Severity-response tables","output/timevarying_pollution_severity_models_with_listing_year","(table|results)[.](html|csv)$")
add("Source concentration maps","output/figures/primary_pollution_study_period_maps","[.](pdf|png)$","Unchanged: pollution source, not lung cohort")
add("Heart/kidney subgroup AJ","output/figures/organ_subgroup_aalen_johansen_cif","[.](pdf|png)$","Unchanged: no lung candidates")
add("SVI validation",file.path(revision_dir,"svi_validation"),"svi_validation_AB[.](pdf|png)$","Unchanged: geographic validation")
manifest <- rbindlist(entries)
fwrite(manifest,file.path(dest,"publication_manifest.csv"))
escape <- function(x) {x<-gsub("&","&amp;",x,fixed=TRUE);x<-gsub("<","&lt;",x,fixed=TRUE);gsub(">","&gt;",x,fixed=TRUE)}
sections <- lapply(unique(manifest$group),function(g) {
  z <- manifest[group==g]
  links <- paste0("<li><a href='../../../",z$path,"'>",escape(basename(z$path)),"</a></li>")
  paste0("<section><h2>",g,"</h2><p>",escape(unique(z$status)),"</p><ul>",paste(links,collapse="\n"),"</ul></section>")
})
writeLines(c("<!doctype html><html><head><meta charset='utf-8'><title>LAS-era publication figures and tables</title><style>body{font:16px Arial;max-width:1100px;margin:40px auto;padding:0 20px;color:#151515}h1{font-size:28px}h2{font-size:20px}section{border-top:1px solid #bbb;padding:14px 0}li{margin:6px 0;overflow-wrap:anywhere}a{color:#075c9c}p{line-height:1.5}</style></head><body>",
  "<h1>LAS-era publication figures and tables</h1><p>Updated October 15, 2026. Final cohort: <strong>1,007,801 candidate-organ observations</strong>, including <strong>42,874 lung candidates</strong>.</p>",
  paste0("<p>",publication_period_note," National AJ quartile thresholds retain the original reference-cohort definition.</p>"),
  "<p>Subgroup estimates marked * have convergence warnings involving the age-under-12 nuisance term. See the model log before interpretation. Map values are model-implied net risks from cause-specific Cox models, not competing-risk cumulative incidence.</p>",
  "<p>This index identifies the current publication displays. Earlier full-era lung proxy analyses, reconstructed LAS/CAS experiments, and obsolete waitlist-average/KM figures remain historical and are not included here. The 2015-2019 three-year competing-risk sensitivity is unchanged because its follow-up ends before the LAS transition.</p>",
  unlist(sections),"</body></html>"),file.path(dest,"publication_index.html"))
writeLines(c("PASS: final cohort and outcome components reconcile.",
  "PASS: all 24 primary baseline/time-varying estimates match their diagnostic refits within 1e-7.",
  "PASS: baseline risk-table counts reconcile to exposure-eligible cohort counts.",
  "PASS: all matched nonlung Table 1 cells equal the pre-update committed values.",
  "PASS: lung map models contain 40,721 candidates and 5,206 adverse events.",
  "PASS: candidate ABC map includes 1,007,801 candidate-organ observations.",
  paste("PASS:",nrow(manifest),"publication artifacts indexed; all files exist and are nonempty.")),file.path(dest,"publication_validation.txt"))
log_revision("Publication checks passed; index generated")
