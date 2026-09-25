source("code/revision_2026/common.R")
env <- load_script_prefix("code/79_timevarying_pollution_cox_svi.R","all_results")
# The primary monthly intervals include the terminal day (tstop = date difference + 1).
env$analysis_dat <- env$analysis_dat[WL_ORG=="KI" & as.numeric(observed_end_date-index_date)+1>3*365.25]
dest <- file.path(revision_dir,"kidney_late_band_audit");dir.create(dest,showWarnings=FALSE)
reference <- read_csv(file.path(revision_dir,"timevarying_diagnostics","followup_band_effects.csv"),show_col_types=FALSE)
warnings_seen <- list();results <- list()
for(i in seq_len(nrow(env$pollutant_specs))) {
  spec <- env$pollutant_specs[i,];log_revision("Kidney late-band warning audit: ",spec$pollutant)
  x <- env$make_month_intervals(env$analysis_dat,spec$exposure_end_date)
  x <- x[tstop>3*365.25]
  x <- env$add_time_updated_scores(x,env$score_sources)
  if(spec$pollutant=="pm25") {
    x <- env$exposure_tables$pm25[x,on=c("zip"="candidate_zip","year","month")]
    x[,pm25_interval_5ug:=pm25_interval_ug_m3/5]
  } else if(spec$pollutant=="o3") {
    x <- env$exposure_tables$o3[x,on=c("zip"="candidate_zip","year","month")]
    x[,o3_interval_10ppb:=o3_interval_ppb/10]
  } else {
    x <- env$exposure_tables$no2_monthly[x,on=c("zip"="candidate_zip","year","month")]
    x <- env$exposure_tables$no2_annual[x,on=c("zip","year")]
    x[,no2_interval_10ppb:=fifelse(!is.na(no2_interval_ppb),no2_interval_ppb,no2_annual_ppb)/10]
  }
  adj <- env$primary_terms("KI")
  needed <- c("tstart","tstop","tv_adverse_event",spec$term,adj$vars)
  x <- x[complete.cases(x[,..needed]) & tstop>tstart]
  form <- as.formula(paste("Surv(tstart,tstop,tv_adverse_event) ~",paste(c(spec$term,adj$terms),collapse=" + ")))
  for(j in 1:2) {
    bounds <- c(3*365.25,5*365.25,Inf)
    band_name <- c("3-5 years","5+ years")[j]
    b <- copy(x[tstop>bounds[j] & tstart<bounds[j+1]])
    b[,tv_adverse_event:=as.integer(tv_adverse_event==1 & tstop<=bounds[j+1])]
    b[,`:=`(tstart=pmax(tstart,bounds[j]),tstop=pmin(tstop,bounds[j+1]))]
    fit <- withCallingHandlers(coxph(form,data=b,ties="efron"),warning=function(w) {
      warnings_seen[[length(warnings_seen)+1L]] <<- data.frame(organ="KI",pollutant=spec$pollutant,phase="followup_band",band=band_name,message=conditionMessage(w),call=paste(deparse(conditionCall(w)),collapse=" "))
      invokeRestart("muffleWarning")
    })
    rr <- broom::tidy(fit,exponentiate=TRUE,conf.int=TRUE) %>% filter(term==spec$term) %>% mutate(organ="KI",pollutant=spec$pollutant,band=band_name,n=uniqueN(b$waitlist_row_id),events=sum(b$tv_adverse_event))
    ref <- reference %>% filter(organ=="KI",pollutant==spec$pollutant,band==band_name)
    stopifnot(nrow(ref)==1,abs(rr$estimate-ref$estimate)<1e-7,rr$n==ref$n,rr$events==ref$events)
    results[[length(results)+1L]] <- rr
    rm(b,fit);gc()
  }
  rm(x);gc()
}
write_csv(bind_rows(warnings_seen),file.path(dest,"captured_warnings.csv"))
write_csv(bind_rows(results),file.path(dest,"reproduced_band_effects.csv"))
