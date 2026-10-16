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
summary <- d[,.(candidates=.N,adverse_events=sum(adverse_event),transplant=sum(component=="transplant"),
  improvement=sum(component=="improvement"),median_waitlist_days=median(followup_days),
  q25=quantile(followup_days,.25),q75=quantile(followup_days,.75)),by=WL_ORG]
summary[, organ_label := organ_names[WL_ORG]]
fwrite(summary,"output/figures/cohort_flow_diagram/cohort_flow_organ_outcomes.csv")

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
fwrite(missing,file.path(out,"cox_model_variables_missingness_long.csv"))
fwrite(missing,file.path(out,"cox_model_variables_missingness_supplement_table.csv"))
gtsave(gt(as.data.frame(missing),groupname_col="Organ") |> fmt_number(columns="Missing (%)",decimals=1) |>
  tab_header(title="Model Covariates and Baseline Missingness") |> tab_source_note(publication_period_note) |>
  tab_source_note("Listing year is not included in either primary model. Baseline models omit organ severity. Missingness here precedes imputation and complete-case model eligibility where raw fields are available; constructed variables may have no missing values by design.") |>
  tab_options(table.font.names="Arial",table.font.size=px(13)),file.path(out,"cox_model_variables_missingness_supplement_table.html"))
writeLines(publication_period_note,file.path(out,"cox_model_variables_missingness_table_note.txt"))

# A flow diagram with the LAS-era restriction explicitly separated from duplicates.
out <- "output/figures/cohort_flow_diagram"
a <- fread(file.path(out,"cohort_flow_accounting.csv"))
collapse <- fread("output/primary_waitlist_period_pollution_cox/overlapping_multilist_collapse_summary.csv")
fmt <- function(x) format(x,big.mark=",",scientific=FALSE,trim=TRUE)
n0 <- a[step=="Eligible registrations before ZCTA exclusion",n]
nm <- a[step=="Excluded: missing residential ZCTA",n]
nz <- a[step=="Included: residential ZCTA available",n]
nd <- sum(collapse$removed_duplicate_registration_rows)
stopifnot(n0-nm==nz,nz-nd==1017540,1017540-9739==nrow(d))
main <- data.table(x=c(5,5,5,5),y=c(17,13.5,10,6.5),w=6.8,h=1.7,
  label=c(paste0("Eligible waitlist registrations, 2005-2025\n",fmt(n0)),
    paste0("Registrations with residential ZCTA\n",fmt(nz)),
    "Deduplicated candidate-organ observations\n1,017,540", "Final analytic cohort\n1,007,801 candidate-organ observations"))
exc <- data.table(x=11.5,y=c(15.3,11.8,8.3),w=4.8,h=c(1.8,2.6,1.8),
  label=c(paste0("Excluded: missing ZCTA\n",fmt(nm)),
    paste0("Duplicate registrations removed: ",fmt(nd),"\n",paste(paste0(collapse$organ_label,": ",fmt(collapse$removed_duplicate_registration_rows)),collapse="\n")),
    "Excluded: lung listing on or after\nMarch 9, 2023\n9,739 candidates"))
summary <- summary[match(c("HR","KI","LI","LU"),WL_ORG)]
org <- data.table(x=c(1.4,5.1,8.8,12.5),y=3.8,w=3.25,h=1.3,
  label=paste0(summary$organ_label,"\n",fmt(summary$candidates)," candidates"))
ends <- data.table(x=org$x,y=1.4,w=3.25,h=2,
  label=paste0("Adverse events: ",fmt(summary$adverse_events),"\nTransplants: ",fmt(summary$transplant),
    "\nMedian waitlist: ",fmt(summary$median_waitlist_days)," days"))
boxes <- rbindlist(list(main,exc,org,ends))
segments <- rbindlist(list(data.table(x=5,xend=5,y=c(16.15,12.65,9.15),yend=c(14.35,10.85,7.35)),
  data.table(x=5,xend=9.1,y=exc$y,yend=exc$y),data.table(x=5,xend=5,y=5.65,yend=5),
  data.table(x=min(org$x),xend=max(org$x),y=5,yend=5),
  data.table(x=org$x,xend=org$x,y=5,yend=4.45),data.table(x=org$x,xend=org$x,y=3.15,yend=2.4)))
p <- ggplot()+geom_segment(data=segments,aes(x=x,xend=xend,y=y,yend=yend),linewidth=.7)+
  geom_rect(data=boxes,aes(xmin=x-w/2,xmax=x+w/2,ymin=y-h/2,ymax=y+h/2),fill="white",color="black",linewidth=.7)+
  geom_text(data=boxes,aes(x,y,label=label),size=5.2,lineheight=1.1)+
  coord_cartesian(xlim=c(-.4,14.4),ylim=c(0,18),expand=FALSE,clip="off")+theme_void()+theme(plot.margin=margin(16,16,16,16))
ggsave(file.path(out,"cohort_flow_diagram.png"),p,width=16,height=18,dpi=320,bg="white")
ggsave(file.path(out,"cohort_flow_diagram.pdf"),p,width=16,height=18,bg="white")
fwrite(boxes,file.path(out,"cohort_flow_diagram_labels.csv"))
writeLines(c(publication_period_note,"The flow applies the lung listing-date restriction after the original candidate-level deduplication. Retained lung follow-up is censored at the transition; 1,274 candidates are affected."),file.path(out,"cohort_flow_note.txt"))
log_revision("Updated cohort flow, outcome components, and covariate missingness")
