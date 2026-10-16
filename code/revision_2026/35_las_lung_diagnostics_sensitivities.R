source("code/revision_2026/common.R")
source("code/publication_cohort.R")
setDTthreads(2)
b <- readRDS(file.path(revision_dir,"cache","las_era_lung","prepared.rds"))$cohort
tv <- readRDS(file.path(revision_dir,"cache","las_era_lung","monthly_intervals.rds"))
update_lung <- function(path,x) {
  old <- fread(path)
  fwrite(rbindlist(list(old[organ!="LU"],as.data.table(x)),fill=TRUE),path)
}
tests <- list(); reproductions <- list(); bands <- list(); eras <- list(); interactions <- list()
warnings <- list()
fit_safe <- function(f,d,context,x=FALSE) {
  withCallingHandlers(coxph(f,data=d,ties="efron",x=x,y=x),warning=function(w) {
    warnings[[length(warnings)+1L]] <<- data.table(context=context,warning=conditionMessage(w))
    invokeRestart("muffleWarning")
  })
}
for(m in c("baseline","tv")) {
  tests <- list(); reproductions <- list(); bands <- list()
  directory <- file.path(revision_dir,if(m=="baseline") "baseline_diagnostics" else "timevarying_diagnostics")
  for(p in names(pollutants)) {
    term <- if(m=="baseline") pollutants[[p]] else c(pm25="pm25_interval_5ug",no2="no2_interval_10ppb",o3="o3_interval_10ppb")[[p]]
    exposure_term <- term
    rhs <- c(term,baseline_adjustment,if(m=="tv") c("las_per10","las_unavailable","las_under12"))
    event <- if(m=="baseline") "adverse_event" else "tv_adverse_event"
    d <- if(m=="baseline") as.data.table(eligible(b,"LU",term)) else copy(tv)
    needed <- c(term,"age","sex","race","zcta_svi_proxy","listing_center",event,if(m=="tv") c("las_per10","las_unavailable","las_under12","tstart","tstop"))
    d <- d[complete.cases(d[,..needed])]
    lhs <- if(m=="baseline") "Surv(followup_days,adverse_event)" else "Surv(tstart,tstop,tv_adverse_event)"
    f <- as.formula(paste(lhs,"~",paste(rhs,collapse="+")))
    log_revision("LAS diagnostics ",m," ",p)
    fit <- fit_safe(f,d,paste(m,p,"primary"),x=TRUE)
    zz <- cox.zph(fit,transform="km",terms=TRUE)
    tests[[p]] <- data.table(organ="LU",pollutant=p,term=rownames(zz$table),as.data.frame(zz$table))
    saveRDS(zz,file.path(revision_dir,"cache",paste0(m,"_zph_LU_",p,".rds")))
    rr <- as.data.table(broom::tidy(fit,exponentiate=TRUE,conf.int=TRUE))[term==exposure_term]
    rr[, `:=`(organ="LU",pollutant=p,n=uniqueN(d$waitlist_row_id),events=sum(d[[event]]),intervals=nrow(d))]
    if(m=="baseline") reproductions[[p]] <- effect_row(fit,term,"baseline_reference_reproduction","LU",p,d) else reproductions[[p]] <- rr
    bounds <- c(0,365.25,3*365.25,5*365.25,Inf)
    for(k in 1:4) {
      if(m=="baseline") {
        dd <- copy(d[followup_days>bounds[k]])
        dd[, `:=`(start=bounds[k],stop=pmin(followup_days,bounds[k+1]),band_event=as.integer(adverse_event==1 & followup_days<=bounds[k+1]))]
        fb <- as.formula(paste("Surv(start,stop,band_event)~",paste(rhs,collapse="+")))
        e <- sum(dd$band_event)
      } else {
        dd <- copy(d[tstop>bounds[k] & tstart<bounds[k+1]])
        dd[,tv_adverse_event:=as.integer(tv_adverse_event==1 & tstop<=bounds[k+1])]
        dd[,`:=`(tstart=pmax(tstart,bounds[k]),tstop=pmin(tstop,bounds[k+1]))]
        fb <- f; e <- sum(dd$tv_adverse_event)
      }
      if(e < if(m=="baseline") 20 else 50) next
      ff <- fit_safe(fb,dd,paste(m,p,"band",k))
      if(m=="baseline") {
        row <- effect_row(ff,term,paste0("baseline_followup_band_",k),"LU",p,dd)
        row$events <- e; row$lower_year <- bounds[k]/365.25; row$upper_year <- bounds[k+1]/365.25
      } else {
        row <- as.data.table(broom::tidy(ff,exponentiate=TRUE,conf.int=TRUE))[term==exposure_term]
        row[,`:=`(organ="LU",pollutant=p,band=c("0-1 years","1-3 years","3-5 years","5+ years")[k],n=uniqueN(dd$waitlist_row_id),events=e)]
      }
      bands[[paste(p,k)]] <- row
    }
    if(m=="baseline") {
      d[, era := factor(ifelse(index_date<as.Date("2017-11-24"),"Era 1","Era 2"))]
      for(e in levels(d$era)) {
        dd <- droplevels(d[era==e]); ff <- fit_safe(f,dd,paste(p,e))
        row <- effect_row(ff,term,paste0("listing_",e),"LU",p,dd); row$adult_heart_only <- FALSE
        eras[[paste(p,e)]] <- row
      }
      adj <- c("age","sex","race","zcta_svi_proxy","strata(listing_center,era)")
      reduced <- fit_safe(as.formula(paste(lhs,"~",paste(c(term,adj),collapse="+"))),d,paste(p,"era reduced"))
      full <- fit_safe(as.formula(paste(lhs,"~",paste(c(paste0(term,":era"),adj),collapse="+"))),d,paste(p,"era full"))
      df <- sum(!is.na(coef(full)))-sum(!is.na(coef(reduced))); stat <- 2*(full$loglik[2]-reduced$loglik[2])
      interactions[[p]] <- data.table(organ="LU",pollutant=p,chisq=stat,df=df,p_value=pchisq(stat,df,lower.tail=FALSE),adult_heart_only=FALSE)
    }
  }
  update_lung(file.path(directory,"schoenfeld_tests.csv"),rbindlist(tests,fill=TRUE))
  update_lung(file.path(directory,"followup_band_effects.csv"),rbindlist(bands,fill=TRUE))
  update_lung(file.path(directory,if(m=="baseline") "baseline_reproduction.csv" else "timevarying_reproduction.csv"),rbindlist(reproductions,fill=TRUE))
}
update_lung(file.path(revision_dir,"baseline_diagnostics","listing_era_effects.csv"),rbindlist(eras))
update_lung(file.path(revision_dir,"baseline_diagnostics","pollutant_by_era_tests.csv"),rbindlist(interactions))
cuts <- fread(file.path(revision_dir,"baseline_diagnostics","policy_era_cutpoints.csv"))
fwrite(cuts[!(organ=="LU" & cut_date=="2023-03-09")],file.path(revision_dir,"baseline_diagnostics","policy_era_cutpoints.csv"))

# Relisting spells and concurrency preserve their distinct sensitivity definitions.
spells <- readRDS(file.path(revision_dir,"cache","spells_prepared.rds"))
spells <- as.data.table(spells)[WL_ORG=="LU" & index_date<=as.Date("2023-03-08")]
spells[observed_end_date>as.Date("2023-03-08"),adverse_event:=0L]
spells[,followup_days:=as.numeric(pmin(observed_end_date,as.Date("2023-03-08"))-index_date)]
ints <- readRDS(file.path(revision_dir,"cache","concurrency.rds"))
regs <- list()
for(p in names(pollutants)) {
  term <- pollutants[[p]]
  f <- as.formula(paste("Surv(followup_days,adverse_event)~",paste(c(term,baseline_adjustment),collapse="+")))
  dd <- eligible(spells,"LU",term)
  ff <- fit_safe(f,dd,paste(p,"spells"))
  rr <- effect_row(ff,term,"separate_spells_model_based","LU",p,dd)
  ff <- coxph(f,data=dd,ties="efron",cluster=PERS_ID)
  rr <- rbind(rr,effect_row(ff,term,"separate_spells_person_clustered","LU",p,dd))
  dd <- eligible(b,"LU",term)
  x <- merge(ints,as.data.table(dd),by.x="primary_id",by.y="waitlist_row_id")
  x <- x[tstart<followup_days]
  x[,tstop:=pmin(tstop,followup_days)]
  x[,interval_event:=as.integer(tstop==followup_days & adverse_event==1)]
  stopifnot(sum(x$interval_event)==sum(dd$adverse_event))
  ff <- fit_safe(as.formula(paste("Surv(tstart,tstop,interval_event)~",paste(c(term,baseline_adjustment,"concurrent_listings","no_open_listing"),collapse="+"))),x,paste(p,"concurrency"))
  row <- effect_row(ff,term,"updated_concurrent_listing_count","LU",p,dd); row$intervals<-nrow(x)
  regs[[p]] <- rbindlist(list(rr,row),fill=TRUE)
  fwrite(regs[[p]],file.path(revision_dir,"registration_models",paste0("LU_",p,".csv")))
  fwrite(as.data.table(broom::tidy(ff,exponentiate=TRUE,conf.int=TRUE)),file.path(revision_dir,"registration_models",paste0("LU_",p,"_concurrency_coefficients.csv")))
}
update_lung(file.path(revision_dir,"registration_models","registration_sensitivity_results.csv"),rbindlist(regs))
old_summary <- fread(file.path(revision_dir,"registration_spell_summary.csv"))
registrations <- as.data.table(readRDS(file.path(revision_dir,"cache","registrations.rds")))
registrations <- registrations[WL_ORG=="LU" & index_date<=as.Date("2023-03-08")]
per_person <- spells[, .(k=.N), by=PERS_ID]
row <- data.table(WL_ORG="LU",people=nrow(b),registrations=nrow(registrations),spells=nrow(spells),
  people_with_multiple_spells=sum(per_person$k>1),extra_spells=sum(per_person$k-1),
  people_with_multiple_registrations=registrations[, .N,by=PERS_ID][N>1,.N])
fwrite(rbindlist(list(old_summary[WL_ORG!="LU"],row)),file.path(revision_dir,"las_era_lung","publication_registration_spell_summary.csv"))
fwrite(rbindlist(warnings,fill=TRUE),file.path(revision_dir,"las_era_lung","updated_diagnostic_warnings.csv"))
log_revision("LAS lung diagnostics, eras, concurrency and relisting sensitivities complete")
