source("code/revision_2026/common.R")
suppressPackageStartupMessages(library(tidyr))
dest <- file.path(revision_dir,"supplement");dir.create(dest,showWarnings=FALSE)
rd <- function(path) read_csv(file.path(revision_dir,path),show_col_types=FALSE)
hr <- function(x,l,u) sprintf("%.2f (%.2f-%.2f)",x,l,u)
p_display <- function(x) ifelse(is.na(x),NA_character_,ifelse(x<.001,"<.001",sprintf("%.3f",x)))
escape <- function(x) {
  x <- gsub("&","&amp;",as.character(x),fixed=TRUE)
  x <- gsub("<","&lt;",x,fixed=TRUE);gsub(">","&gt;",x,fixed=TRUE)
}
export <- function(d,name,caption,note,era_context=NULL) {
  if(any(grepl("[a]",as.character(unlist(d)),fixed=TRUE))) note <- paste(note,"[a] A sparse categorical nuisance coefficient triggered a possible infinite-coefficient warning. The pollutant coefficient remained finite; interpret the flagged model cautiously. No race categories were silently combined to remove warnings.")
  note <- paste(note,"Abbreviations: ACS, American Community Survey; CDC, Centers for Disease Control and Prevention; CI, confidence interval; HR, hazard ratio; SE, standard error; SVI, Social Vulnerability Index; ZCTA, ZIP Code Tabulation Area; PM2.5, fine particulate matter; NO2, nitrogen dioxide; O3, ozone. Organ codes: HR, heart; KI, kidney; LI, liver; LU, lung.")
  write_csv(d,file.path(dest,paste0(name,".csv")))
  write_tsv(d,file.path(dest,paste0(name,".tsv")))
  header <- paste0("<tr>",paste0("<th>",escape(names(d)),"</th>",collapse=""),"</tr>")
  rows <- apply(as.data.frame(d),1,function(x) paste0("<tr>",paste0("<td>",escape(x),"</td>",collapse=""),"</tr>"))
  context_html <- ""
  if(!is.null(era_context)) {
    context_header <- paste0("<tr>",paste0("<th>",escape(names(era_context)),"</th>",collapse=""),"</tr>")
    context_rows <- apply(as.data.frame(era_context),1,function(x) paste0("<tr>",paste0("<td>",escape(x),"</td>",collapse=""),"</tr>"))
    context_html <- paste0("<h3>Policy-era definitions and rationale</h3><table>",context_header,paste(context_rows,collapse=""),"</table>")
  }
  writeLines(c("<!doctype html><html><head><meta charset='utf-8'><style>body{font:15px Arial;margin:32px;color:#111}table{border-collapse:collapse}th,td{padding:8px;border-bottom:1px solid #bbb;text-align:left;overflow-wrap:anywhere}th{border-top:2px solid #111;border-bottom:2px solid #111}p{max-width:1200px}</style></head><body>",paste0("<h2>",escape(caption),"</h2><table>",header,paste(rows,collapse=""),"</table><p>",escape(note),"</p>",context_html,"</body></html>")),file.path(dest,paste0(name,".html")))
}
ref <- rd("baseline_diagnostics/baseline_reproduction.csv") %>% mutate(analysis="Primary baseline")
reg <- rd("registration_models/registration_sensitivity_results.csv")
recent_folder <- if(file.exists(file.path(revision_dir,"finegray_2015_2019","validation.txt"))) "finegray_2015_2019" else "recent_finegray"
recent_window <- if(recent_folder=="finegray_2015_2019") "2015-2019" else "2020-2022"
fg <- rd(paste0(recent_folder,"/finegray_results.csv"))
cs <- rd(paste0(recent_folder,"/matched_recent_cause_specific.csv"))
recent_warnings <- bind_rows(fg,cs)
if("convergence_note" %in% names(recent_warnings)) {
  recent_warnings <- recent_warnings %>% filter(!is.na(convergence_note),convergence_note!="") %>%
    transmute(analysis,organ,pollutant,message=convergence_note)
} else recent_warnings <- tibble(analysis=character(),organ=character(),pollutant=character(),message=character())
fg <- fg %>% select(-any_of("convergence_note"))
cs <- cs %>% select(-any_of("convergence_note"))
audit_base <- rd("baseline_warning_audit/captured_warnings.csv")
audit_tv <- bind_rows(rd("thoracic_warning_audit/captured_warnings.csv"),rd("kidney_late_band_audit/captured_warnings.csv"))
warning_summary <- function(w,keys) w %>% group_by(across(all_of(keys))) %>% summarise(convergence_note=paste(unique(message),collapse=" | "),.groups="drop")
all <- bind_rows(ref,reg,fg,cs) %>% mutate(estimate_ci=hr(hazard_ratio,conf_low,conf_high))
stopifnot(nrow(ref)==12,nrow(reg)==36,nrow(fg)==12,nrow(cs)==12,all(is.finite(all$hazard_ratio)))
w <- audit_base %>% mutate(analysis=recode(analysis,primary_baseline="Primary baseline",recent_FineGray="recent_2020_2022_3yr_FineGray",recent_cause_specific="recent_2020_2022_3yr_cause_specific"))
w <- bind_rows(w,recent_warnings)
all <- left_join(all,warning_summary(w,c("analysis","organ","pollutant")),by=c("analysis","organ","pollutant")) %>% mutate(convergence_note=coalesce(convergence_note,""),estimate_ci=if_else(convergence_note!="",paste0(estimate_ci," [a]"),estimate_ci))
all <- all %>% mutate(Organ=recode(organ,!!!organ_names),Pollutant=recode(pollutant,pm25="PM2.5",no2="NO2",o3="O3"))
all <- all %>% mutate(analysis=recode(analysis,
  separate_spells_model_based="Separate spells",
  separate_spells_person_clustered="Separate spells (clustered SE)",
  updated_concurrent_listing_count="Concurrent listings",
  recent_2020_2022_3yr_FineGray="Recent cohort: Fine-Gray",
  recent_2020_2022_3yr_cause_specific="Recent cohort: cause-specific",
  recent_2015_2019_3yr_FineGray="Recent cohort: Fine-Gray",
  recent_2015_2019_3yr_cause_specific="Recent cohort: cause-specific"))
export(all %>% mutate(p_value=p_display(p_value)) %>% select(Organ,Pollutant,analysis,n,people,events,estimate_ci,p_value,convergence_note),"sensitivity_models_long",
"Supplemental Table. Registration and recent-cohort sensitivity analyses",
paste0("Estimates are HRs (95% CIs), except Fine-Gray estimates, which are subdistribution HRs. PM2.5: per 5 ug/m3; NO2 and O3: per 10 ppb. All use preceding-year exposure, age, sex, race, ACS vulnerability proxy and center strata. Concurrency adds updated listing count and a no-open-listing indicator. Separate spells reset baseline after gaps. Recent-cohort models restrict listings to ",recent_window," with follow-up capped at 3 years; they do not use the full-period cohort. Pandemic-era follow-up is retained. Person-clustered variance is reported separately for spells and used for Fine-Gray pseudo-records."))
wide <- all %>% select(Organ,Pollutant,analysis,estimate_ci) %>% pivot_wider(names_from=analysis,values_from=estimate_ci)
export(wide,"sensitivity_models_wide","Supplemental Table. Sensitivity effect estimates",
paste0("HR (95% CI); Fine-Gray column reports subdistribution HR. Recent-cohort columns use ",recent_window," listings with 3-year follow-up, including pandemic-era follow-up, and have different eligibility than full-period columns. See accompanying long-format table for sample sizes, events, adjustment and units."))
if(recent_folder=="finegray_2015_2019") {
  export(all %>% filter(grepl("^Recent cohort:",analysis)) %>% select(Organ,Pollutant,analysis,n,events,estimate_ci,p_value,convergence_note) %>% mutate(p_value=p_display(p_value)),
    "competing_risks_2015_2019","Supplemental Table. Three-year outcomes among candidates listed in 2015-2019",
    "All eligible listings from January 1, 2015 through December 31, 2019, using the primary collapsed candidate-organ cohort; no subsampling. Exposure is the preceding-year average. Adjustment: age, sex, race, ACS vulnerability proxy and center strata; no listing-year or organ-score term. Follow-up is capped at 3*365.25 days and can extend into the pandemic. Fine-Gray models treat transplant/improvement as competing events and other exits as censoring; center-specific censoring weights and person-clustered sandwich SEs account for pseudo-record expansion. Matched cause-specific models censor competing events. Estimates are subdistribution HRs for Fine-Gray and cause-specific HRs otherwise, with 95% CIs. PM2.5: per 5 ug/m3; NO2/O3: per 10 ppb. These are distinct estimands and neither HR is a risk ratio.")
}
eras <- rd("baseline_diagnostics/listing_era_effects.csv")
cuts <- rd("baseline_diagnostics/policy_era_cutpoints.csv")
era_context <- tibble::tribble(
  ~Organ, ~`New era begins`, ~`Policy change and rationale`, ~Source,
  "Heart", "2018-10-18", "Revised adult heart allocation with six urgency statuses changed medical priority and organ sharing. Heart era analyses include adults only.", "https://hrsa.unos.org/media/opqmziro/data_report_heart_committee_5yr_rpt2_042024.pdf",
  "Kidney", "2014-12-04", "Kidney Allocation System introduced changes to matching and priority, including credit for prelisting dialysis time, potentially changing access and waitlist outcomes.", "https://hrsa.unos.org/media/1235/kas_faqs.pdf",
  "Kidney", "2021-03-15", "Circle-based kidney allocation replaced donation-service-area boundaries, changing geographic access to deceased-donor transplantation.", "https://pubmed.ncbi.nlm.nih.gov/37196709/",
  "Liver", "2016-01-11", "MELD-Na incorporated serum sodium into medical-priority scoring, changing candidate ranking.", "https://unos.org/news/policy-and-system-changes-effective-january-11-2016-adding-serum-sodium-to-meld-calculation/",
  "Liver", "2020-02-04", "Acuity-circle allocation replaced donation-service-area and regional boundaries with distance-based distribution, changing geographic access.", "https://unos.org/news/pre-imp-notice-liver-intestinal-dist-acuity-circles-feb-4-2020/",
  "Liver", "2023-07-13", "MELD 3.0 revised medical-priority scoring, including sex and albumin, changing candidate ranking and addressing sex disparities.", "https://hrsa.unos.org/media/fyxhlkp5/improving-liver-allocation-meld-30-faq.pdf",
  "Lung", "2017-11-24", "Initial geographic distribution changed from donation service areas to a 250-nautical-mile radius, broadening access.", "https://hrsa.unos.org/media/jlijyukx/thoracic_publiccomment_distribution_20180122.pdf",
  "Lung", "2023-03-09", "Continuous distribution using the composite allocation score (CAS) replaced LAS-based allocation, changing the combination of medical urgency, benefit and access factors.", "https://hrsa.unos.org/news/early-lung-monitoring-report-shows-increase-in-transplants-following-implementation-of-lung-continuous-distribution/"
)
stopifnot(identical(sort(as.character(cuts$cut_date)),sort(era_context$`New era begins`)))
era_note <- paste("Era boundaries represent selected major changes in allocation priority or geographic access, not every policy revision.",
  "The first era starts at study entry in 2005; each subsequent era starts on the date below and ends before the next boundary, or at the pollutant-specific study end.",
  "Era is assigned at listing; follow-up is not split or censored at a policy transition. These analyses do not estimate causal policy effects.",
  "Adjustment: age, sex, race, ACS vulnerability proxy and transplant-center strata; no listing-year or organ-score term.",
  "MELD, Model for End-Stage Liver Disease; MELD-Na, MELD incorporating sodium; LAS, lung allocation score; CAS, composite allocation score.")
export(era_context,"policy_era_rationale","Supplemental Table. Policy-era definitions and rationale",era_note)
era_labels <- bind_rows(lapply(names(organ_names),function(o) {
  x <- as.character(cuts$cut_date[cuts$organ==o]);data.frame(organ=o,analysis=paste0("listing_Era ",seq_len(length(x)+1)),listing_start=c("Study start",x),listing_end_exclusive=c(x,"Study end"))
}))
eras <- left_join(eras,era_labels,by=c("organ","analysis")) %>% mutate(estimate_ci=hr(hazard_ratio,conf_low,conf_high))
eras <- left_join(eras,warning_summary(audit_base,c("analysis","organ","pollutant")),by=c("analysis","organ","pollutant")) %>% mutate(convergence_note=coalesce(convergence_note,""),estimate_ci=if_else(convergence_note!="",paste0(estimate_ci," [a]"),estimate_ci))
export(eras %>% mutate(p_value=p_display(p_value)) %>% select(organ,pollutant,listing_start,listing_end_exclusive,n,events,estimate_ci,p_value,adult_heart_only,convergence_note),"policy_era_models","Supplemental Table. Listing-era-specific baseline Cox models",
paste("Intervals include their starting date and exclude the next era start. PM2.5 HRs are per 5 ug/m3; NO2 and O3 HRs are per 10 ppb. Separate-era point estimates allow nuisance coefficients to differ.",era_note),era_context=era_context)
export(rd("baseline_diagnostics/pollutant_by_era_tests.csv") %>% mutate(chisq=round(chisq,2),p_value=p_display(p_value)),"policy_era_interactions","Supplemental Table. Pollutant-by-listing-era interaction tests",
paste("Likelihood-ratio tests compare a common pollutant coefficient with era-specific pollutant coefficients, using center-by-era strata and common nuisance coefficients. The null is equal pollutant associations across all listing eras; these are omnibus tests, not pairwise comparisons. P values are not adjusted for the 12 exploratory interaction tests. Read alongside era-specific HRs and confidence intervals.",era_note),era_context=era_context)
ph <- bind_rows(rd("baseline_diagnostics/schoenfeld_tests.csv") %>% mutate(model="Baseline"),rd("timevarying_diagnostics/schoenfeld_tests.csv") %>% mutate(model="Time-varying"))
stopifnot(nrow(distinct(ph,model,organ,pollutant))==24)
export(ph %>% mutate(chisq=round(chisq,2),p=p_display(p)),"proportional_hazards_tests","Supplemental Table. Proportional-hazards diagnostics","Scaled Schoenfeld tests use the Kaplan-Meier time transform. Tests are reported for all nonstratified terms and globally. Small P values indicate evidence against a constant coefficient over time; they are not evidence that an exposure association is absent. Time-band estimates and residual plots accompany these tests.")
bb <- rd("baseline_diagnostics/followup_band_effects.csv") %>% mutate(analysis=sub("baseline_followup_band_","baseline_band_",analysis)) %>% left_join(warning_summary(audit_base,c("analysis","organ","pollutant")),by=c("analysis","organ","pollutant")) %>% mutate(convergence_note=coalesce(convergence_note,""),estimate_ci=paste0(hr(hazard_ratio,conf_low,conf_high),if_else(convergence_note!=""," [a]","")),p_value=p_display(p_value))
export(bb %>% select(organ,pollutant,lower_year,upper_year,n,events,estimate_ci,p_value,convergence_note),"baseline_followup_bands","Supplemental Table. Baseline exposure associations by time since listing","Bands: 0-1, 1-3, 3-5 and 5+ years. Candidate counts are those entering each band. Conditional risk sets change over follow-up; these are descriptive diagnostics, not causal comparisons between bands.")
tb <- rd("timevarying_diagnostics/followup_band_effects.csv") %>% left_join(warning_summary(audit_tv,c("organ","pollutant","band")),by=c("organ","pollutant","band")) %>% mutate(convergence_note=coalesce(convergence_note,""),estimate_ci=paste0(hr(estimate,conf.low,conf.high),if_else(convergence_note!=""," [a]","")),p_value=p_display(p.value))
export(tb %>% select(organ,pollutant,band,n,events,estimate_ci,p_value,convergence_note),"timevarying_followup_bands","Supplemental Table. Time-varying exposure associations by time since listing","Monthly start-stop intervals are clipped at each band boundary without assigning a future event to an earlier band. Primary time-varying covariates and pollutant units are retained. Bands with fewer than 50 adverse events are omitted.")
diagnostic_bands <- bind_rows(
  bb %>% transmute(organ,pollutant,model="Baseline",band=case_when(lower_year==0~"0-1 years",lower_year==1~"1-3 years",lower_year==3~"3-5 years",TRUE~"5+ years"),estimate_ci),
  tb %>% transmute(organ,pollutant,model="Time-varying",band,estimate_ci)
) %>% pivot_wider(names_from=band,values_from=estimate_ci)
diagnostic_summary <- ph %>% filter(grepl("^(pm25|no2|o3)_",term)) %>%
  transmute(organ,pollutant,model,`Pollutant PH P value`=p_display(p)) %>%
  left_join(diagnostic_bands,by=c("organ","pollutant","model")) %>%
  arrange(match(organ,c("HR","KI","LI","LU")),match(model,c("Baseline","Time-varying")),match(pollutant,c("pm25","no2","o3"))) %>%
  mutate(Organ=recode(organ,!!!organ_names),Model=recode(model,Baseline="Prelisting",`Time-varying`="Time-varying"),
    Pollutant=recode(pollutant,pm25="PM2.5",no2="NO2",o3="O3")) %>%
  select(Organ,Model,Pollutant,`Pollutant PH P value`,`0-1 years`,`1-3 years`,`3-5 years`,`5+ years`) %>%
  mutate(across(everything(),~replace_na(as.character(.x),"Not estimated")))
stopifnot(nrow(diagnostic_summary)==24L,!anyDuplicated(diagnostic_summary[c("Organ","Model","Pollutant")]),
  all(ph$p[ph$term=="GLOBAL"]<.001))
export(diagnostic_summary,"ph_diagnostics_summary","Supplemental Table. Proportional-hazards diagnostics and follow-up-specific pollution associations",
  paste("Time-interval columns report adjusted HRs (95% CIs), per 5 ug/m3 higher PM2.5 or 10 ppb higher NO2/O3. Prelisting models use exposure in the preceding year; time-varying models update exposure during follow-up.",
    "PH (proportional hazards) P values are pollutant-specific scaled Schoenfeld residual tests using the Kaplan-Meier time transform; a small P value indicates evidence of a changing coefficient, not the magnitude of that change. P values are not multiplicity-adjusted. Global PH tests were P<.001 in all 24 models, including models with nonsignificant pollutant-specific tests.",
    "Separate models within 0-1, 1-3, 3-5 and 5+ years retain the corresponding primary adjustment set and center strata. Intervals are (0,1], (1,3], (3,5] and (5,infinity) years. Their changing risk sets and imprecision at later times should be considered; these are not causal comparisons between intervals.",
    "Not estimated indicates fewer than 20 adverse events for a baseline interval or fewer than 50 for a time-varying interval. Updating exposure does not itself allow its coefficient to vary over time. Residual plots were also generated with natural cubic spline smooths (4 degrees of freedom) and approximate pointwise 95% bands."))
export(bind_rows(audit_base,audit_tv,recent_warnings),"convergence_warning_audit","Supplemental Table. Model convergence warnings","Warnings are reported explicitly. Baseline sensitivity warnings identify the affected coefficient names. The time-varying audit identifies affected follow-up bands. The updated recent-cohort warnings are included with their analysis identifiers; full coefficient files are retained in the corresponding results directory. These are distinct from Schoenfeld proportional-hazards tests.")
export(rd("svi_validation/validation_results.csv") %>% mutate(across(-c(cdc_measure,n_zctas),~round(.x,3))),"svi_validation","Supplemental Table. ACS vulnerability proxy validation against CDC SVI","One observation per matched ZCTA, observed 2022 ACS vintage. Larger values indicate greater vulnerability. Rank correlations, quartile agreement and linearly weighted kappa assess convergent validity, not interchangeability. Bootstrap intervals resample ZCTAs and do not model spatial dependence. Missing proxy components retain the primary median-imputation method.")
base <- baseline_data()
ints <- readRDS(file.path(revision_dir,"cache","concurrency.rds"))
entry <- ints[,.(entry_count=concurrent_listings[which.min(tstart)],maximum_count=max(concurrent_listings)),by=primary_id]
entry <- merge(entry,base[,.(primary_id=waitlist_row_id,WL_ORG)],by="primary_id")
write_csv(entry[,.N,by=.(WL_ORG,entry_count)][order(WL_ORG,entry_count)],file.path(dest,"entry_listing_count_distribution.csv"))
export(rd("registration_spell_summary.csv"),"registration_spell_counts","Supplemental Table. Registration and separate-spell counts","Counts describe the reconstructed cohort before pollutant availability and complete-case restrictions. Same-day touching or overlapping registrations are combined; a positive gap begins a new spell. Model-specific eligible counts are in the sensitivity table.")
writeLines(c(paste("Primary reference git commit:",system2("git",c("rev-parse","HEAD"),stdout=TRUE)),paste("Generated:",Sys.time()),capture.output(sessionInfo())),file.path(dest,"provenance.txt"))
