source("code/revision_2026/common.R")
dat <- baseline_data()
cohort_start <- getOption("revision.recent_start","2020-01-01")
cohort_end <- getOption("revision.recent_end","2022-12-31")
cohort_tag <- paste0(substr(cohort_start,1,4),"_",substr(cohort_end,1,4))
dest<-file.path(revision_dir,getOption("revision.recent_dir","recent_finegray"));dir.create(dest,showWarnings=FALSE)
stopifnot(max(dat$observed_end_date,na.rm=TRUE)>=as.Date("2025-12-31"))
results<-list()
for(org in names(organ_names)) for(p in names(pollutants)) {
  file<-file.path(dest,paste0(org,"_",p,".csv"))
  if(file.exists(file)) { results[[length(results)+1L]]<-read_csv(file,show_col_types=FALSE);next }
  term<-pollutants[[p]]
  d<-eligible(dat[index_date>=as.Date(cohort_start) & index_date<=as.Date(cohort_end)],org,term)
  stopifnot(nrow(d)>0L,!anyNA(d$transplant_or_improvement),!any(d$adverse_event==1 & d$transplant_or_improvement==1))
  d$ftime<-pmin(d$followup_days,3*365.25)
  d$event<-factor(ifelse(d$followup_days>3*365.25,0,ifelse(d$adverse_event==1,1,ifelse(d$transplant_or_improvement==1,2,0))),levels=0:2,labels=c("censor","adverse","transplant_improvement"))
  log_revision("Fine-Gray ",org," ",p," candidates=",nrow(d)," adverse=",sum(d$event=="adverse"))
  # Center-specific censoring weights AND center-stratified subdistribution hazards.
  centers<-split(d,d$listing_center,drop=TRUE)
  pieces<-lapply(names(centers),function(center) {
    dd<-centers[[center]]
    if(!any(dd$event=="adverse")) return(NULL)
    ff<-finegray(as.formula(paste("Surv(ftime,event) ~",paste(c(term,"age","sex","race","zcta_svi_proxy","PERS_ID"),collapse=" + "))),data=dd,etype="adverse")
    ff$listing_center<-center
    ff
  })
  fg<-as.data.frame(rbindlist(pieces,use.names=TRUE));rm(pieces,centers);gc()
  log_revision("Fine-Gray expanded rows=",nrow(fg))
  warnings <- character()
  fit<-withCallingHandlers(coxph(as.formula(paste("Surv(fgstart,fgstop,fgstatus) ~",paste(c(term,baseline_adjustment),collapse=" + "))),
             data=fg,weights=fgwt,cluster=PERS_ID,ties="breslow",x=FALSE,y=FALSE),warning=function(w) {
               warnings <<- c(warnings,conditionMessage(w))
             })
  rr<-effect_row(fit,term,paste0("recent_",cohort_tag,"_3yr_FineGray"),org,p,d)
  rr$convergence_note <- paste(unique(warnings),collapse=" | ")
  write_csv(broom::tidy(fit),file.path(dest,paste0(org,"_",p,"_coefficients.csv")))
  rr$events<-sum(d$event=="adverse");rr$competing_events<-sum(d$event=="transplant_improvement")
  rr$expanded_rows<-nrow(fg);rr$variance<-"person_clustered_sandwich_for_FineGray_expansion"
  rr$centers_with_events<-length(unique(fg$listing_center))
  write_csv(rr,file);results[[length(results)+1L]]<-rr
  rm(fg,fit);gc()
}
write_csv(bind_rows(results),file.path(dest,"finegray_results.csv"))
writeLines(c(paste0("Listings from ",cohort_start," through ",cohort_end," inclusive; all eligible candidates; follow-up capped at 3*365.25 days."),
 "Listing-era restriction does not exclude pandemic-era follow-up. The collapsed primary candidate-organ cohort is retained.",
 "Prior-year exposure; age, sex, race, ACS proxy; center-specific censoring weights and center-stratified Fine-Gray baseline hazards.",
 "Transplant/improvement competing; other exits censored. No subsampling. Centers without adverse events contribute no coefficient information.",
 "Breslow handling of tied event times; person-clustered sandwich variance accounts for Fine-Gray pseudo-record expansion."),file.path(dest,"methods.txt"))
