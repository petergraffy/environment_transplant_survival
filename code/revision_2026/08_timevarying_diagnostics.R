source("code/revision_2026/common.R")
options(warn=1)
env <- load_script_prefix("code/79_timevarying_pollution_cox_svi.R","all_results")
dest <- file.path(revision_dir,"timevarying_diagnostics");dir.create(dest,showWarnings=FALSE)
reference <- read_csv(file.path(env$out_dir,"timevarying_pollution_cox_svi_results.csv"),show_col_types=FALSE)
tests <- list(); effects <- list(); bands <- list()
for(i in seq_len(nrow(env$pollutant_specs))) {
  spec <- env$pollutant_specs[i,]
  log_revision("Diagnostic intervals: ",spec$pollutant)
  intervals <- env$make_month_intervals(env$analysis_dat,spec$exposure_end_date)
  intervals <- env$add_time_updated_scores(intervals,env$score_sources)
  if(spec$pollutant=="pm25") {
    intervals <- env$exposure_tables$pm25[intervals,on=c("zip"="candidate_zip","year","month")]
    intervals[,pm25_interval_5ug:=pm25_interval_ug_m3/5]
  } else if(spec$pollutant=="o3") {
    intervals <- env$exposure_tables$o3[intervals,on=c("zip"="candidate_zip","year","month")]
    intervals[,o3_interval_10ppb:=o3_interval_ppb/10]
  } else {
    intervals <- env$exposure_tables$no2_monthly[intervals,on=c("zip"="candidate_zip","year","month")]
    intervals <- env$exposure_tables$no2_annual[intervals,on=c("zip","year")]
    intervals[,no2_interval_10ppb:=fifelse(!is.na(no2_interval_ppb),no2_interval_ppb,no2_annual_ppb)/10]
  }
  for(org in names(organ_names)) {
    prefix <- paste0(org,"_",spec$pollutant)
    adj <- env$primary_terms(org)
    needed <- c("tstart","tstop","tv_adverse_event",spec$term,adj$vars)
    d <- intervals[WL_ORG==org & complete.cases(intervals[,..needed]) & tstop>tstart]
    form <- as.formula(paste("Surv(tstart,tstop,tv_adverse_event) ~",paste(c(spec$term,adj$terms),collapse=" + ")))
    log_revision("Time-varying PH diagnostics: ",prefix," rows=",nrow(d))
    fit <- coxph(form,data=d,ties="efron",x=TRUE,y=TRUE)
    rr <- broom::tidy(fit,exponentiate=TRUE,conf.int=TRUE) %>% filter(term==spec$term) %>% mutate(organ=org,pollutant=spec$pollutant,intervals=nrow(d),n=uniqueN(d$waitlist_row_id),events=sum(d$tv_adverse_event))
    ref <- reference %>% filter(organ==.env$org,pollutant==spec$pollutant)
    stopifnot(nrow(ref)==1L,abs(rr$estimate-ref$hazard_ratio)<1e-7,rr$n==ref$candidate_organ_episodes,rr$events==ref$adverse_events)
    effects[[prefix]] <- rr
    z <- cox.zph(fit,transform="km",terms=TRUE)
    zz <- as.data.frame(z$table) %>% rownames_to_column("term") %>% mutate(organ=org,pollutant=spec$pollutant)
    tests[[prefix]] <- zz
    write_csv(zz,file.path(dest,paste0(prefix,"_tests.csv")))
    saveRDS(z,file.path(revision_dir,"cache",paste0("tv_zph_",prefix,".rds")))
    for(ext in c("png","pdf")) {
      path <- file.path(dest,paste0(prefix,"_schoenfeld.",ext))
      if(ext=="png") png(path,width=2200,height=1500,res=220) else pdf(path,width=10,height=6.8)
      plot(z,var=spec$term,ylab="Pollutant log hazard ratio",xlab="Days since listing",main=paste(organ_names[[org]],toupper(spec$pollutant)))
      abline(h=coef(fit)[spec$term],lty=2,col="gray40")
      dev.off()
    }
    rm(fit,z);gc()
    bounds <- c(0,365.25,3*365.25,5*365.25,Inf)
    for(j in 1:4) {
      b <- copy(d[tstop>bounds[j] & tstart<bounds[j+1]])
      b[,tv_adverse_event:=as.integer(tv_adverse_event==1 & tstop<=bounds[j+1])]
      b[,`:=`(tstart=pmax(tstart,bounds[j]),tstop=pmin(tstop,bounds[j+1]))]
      if(sum(b$tv_adverse_event)<50L) next
      log_revision("Time-band fit: ",prefix," band=",j," events=",sum(b$tv_adverse_event))
      f <- coxph(form,data=b,ties="efron")
      bands[[paste(prefix,j)]] <- broom::tidy(f,exponentiate=TRUE,conf.int=TRUE) %>% filter(term==spec$term) %>% mutate(organ=org,pollutant=spec$pollutant,band=c("0-1 years","1-3 years","3-5 years","5+ years")[j],n=uniqueN(b$waitlist_row_id),events=sum(b$tv_adverse_event))
      rm(b,f);gc()
    }
    write_csv(bind_rows(effects),file.path(dest,"timevarying_reproduction.csv"))
    write_csv(bind_rows(tests),file.path(dest,"schoenfeld_tests.csv"))
    write_csv(bind_rows(bands),file.path(dest,"followup_band_effects.csv"))
    rm(d);gc()
  }
  rm(intervals);gc()
}
