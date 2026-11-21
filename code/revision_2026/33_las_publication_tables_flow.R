source("code/revision_2026/common.R")
source("code/publication_cohort.R")
source("code/revision_2026/outcome_components.R")
library(gt)
setDTthreads(2)
dest <- file.path(revision_dir,"las_era_lung")
d <- publication_baseline_data()
d[, lung_las_era := WL_ORG == "LU"]
components <- primary_outcome_components(d)
d[, component := NULL]
d <- merge(d, components[,.(PERS_ID,WL_ORG,component)],by=c("PERS_ID","WL_ORG"),all.x=TRUE)
stopifnot(nrow(d)==1007801L, !anyNA(d$component))
counts <- dcast(d[,.(n=.N),by=.(WL_ORG,component)], WL_ORG~component,value.var="n",fill=0)
counts[, N := rowSums(.SD),.SDcols=setdiff(names(counts),"WL_ORG")]
fwrite(counts,file.path(revision_dir,"death_only_sensitivity","outcome_components_by_organ.csv"))
fwrite(counts,file.path(dest,"publication_outcomes_by_organ.csv"))

# Covariate missingness is measured before model-specific complete-case selection.
vars <- c(age="Age at listing",sex="Sex",race="Race",zcta_svi_proxy="ACS vulnerability proxy",listing_center="Listing center",
  pm25_prior_ug_m3="PM2.5 prior-365-day exposure",no2_prior_ppb="NO2 prior-12-month exposure",o3_prior_ppb="O3 prior-365-day exposure")
rows <- list()
add <- function(org,var,label,role,handling) {
  x <- d[WL_ORG==org][[var]]
  missing <- is.na(x) | as.character(x) %in% c("","Missing")
  rows[[length(rows)+1L]] <<- data.table(Organ=organ_names[[org]],Variable=var,Description=label,Model=role,
    N=length(x),Missing=sum(missing),`Missing (%)`=100*mean(missing),Handling=handling)
}
for(org in names(organ_names)) {
  for(v in names(vars)) add(org,v,vars[[v]],if(grepl("prior",v)) "Baseline exposure" else "Baseline and time-varying adjustment",
    if(v=="listing_center") "Stratified baseline hazard" else if(v %in% c("sex","race")) "Categorical; explicit missing category where present" else
    if(v=="zcta_svi_proxy") "Six component ranks averaged; source components median-imputed" else "Continuous; required for model fitting")
  if(org %in% c("HR","LI")) add(org,"organ_score",if(org=="HR") "US-CRS" else "MELD/PELD", "Time-varying severity (baseline availability shown)","Existing organ-specific score construction and updates retained")
  if(org=="HR") for(v in c("hr_albumin_observed","hr_bilirubin_observed","hr_sodium_observed","hr_bnp_observed","hr_bnp_type","hr_egfr","hr_durable_lvad","hr_short_mcs"))
    add(org,v,v,"US-CRS component at listing",if(grepl("observed",v)) "Observed missingness; laboratory values median-imputed for the score" else "Derived or categorical component; construction retained from primary analysis")
  if(org=="KI") for(v in c("kidney_no_dialysis_time","kidney_dialysis_years","kidney_diabetes")) add(org,v,v,"Time-varying severity (baseline availability shown)","Separate clinical covariate; dialysis duration updated during follow-up")
  if(org=="LU") add(org,"baseline_las","Official LAS available at listing","Time-varying severity (baseline availability shown)",
    "Monthly latest observed LAS/10; carry forward, no future backfill. Missing filled with 37.80857; missing and age<12 indicators included")
}
missing <- rbindlist(rows)
out <- "output/tables/primary_cox_model_variables_missingness"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
fwrite(missing,file.path(out,"cox_model_variables_missingness_long.csv"))
fwrite(missing,file.path(out,"cox_model_variables_missingness_supplement_table.csv"))
gtsave(gt(as.data.frame(missing),groupname_col="Organ") |> fmt_number(columns="Missing (%)",decimals=1) |>
  tab_header(title="Model Covariates and Baseline Missingness") |> tab_source_note(publication_period_note) |>
  tab_source_note("Listing year is not included in either primary model. Baseline models omit organ severity. Missingness here precedes imputation and complete-case model eligibility where raw fields are available; constructed variables may have no missing values by design.") |>
  tab_options(table.font.names="Arial",table.font.size=px(13)),file.path(out,"cox_model_variables_missingness_supplement_table.html"))
writeLines(publication_period_note,file.path(out,"cox_model_variables_missingness_table_note.txt"))

log_revision("Updated outcome components and covariate missingness")
