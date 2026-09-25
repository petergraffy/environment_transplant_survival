source("code/revision_2026/common.R")
base <- baseline_data()
spells <- readRDS(file.path(revision_dir,"cache","spells_prepared.rds"))
ints <- readRDS(file.path(revision_dir,"cache","concurrency.rds"))
dest <- file.path(revision_dir,"registration_models"); dir.create(dest,showWarnings=FALSE)
results <- list()
for(org in names(organ_names)) for(p in names(pollutants)) {
  term <- pollutants[[p]]
  log_revision("Registration sensitivity: ",org," ",p)
  d <- eligible(spells,org,term)
  form <- as.formula(paste("Surv(followup_days,adverse_event) ~",paste(c(term,baseline_adjustment),collapse=" + ")))
  fit <- coxph(form,data=d,ties="efron")
  rr <- effect_row(fit,term,"separate_spells_model_based",org,p,d)
  fit <- coxph(form,data=d,ties="efron",cluster=PERS_ID)
  rr <- bind_rows(rr,effect_row(fit,term,"separate_spells_person_clustered",org,p,d))
  d <- eligible(base,org,term)
  x <- merge(ints,as.data.table(d),by.x="primary_id",by.y="waitlist_row_id",sort=FALSE)
  x[, interval_event := as.integer(tstop==followup_days & adverse_event==1L)]
  stopifnot(sum(x$interval_event)==sum(d$adverse_event),uniqueN(x$primary_id)==nrow(d))
  fit <- coxph(as.formula(paste("Surv(tstart,tstop,interval_event) ~",paste(c(term,baseline_adjustment,"concurrent_listings","no_open_listing"),collapse=" + "))),data=x,ties="efron")
  r <- effect_row(fit,term,"updated_concurrent_listing_count",org,p,d)
  r$intervals <- nrow(x)
  rr <- bind_rows(rr,r)
  write_csv(broom::tidy(fit,exponentiate=TRUE,conf.int=TRUE),file.path(dest,paste0(org,"_",p,"_concurrency_coefficients.csv")))
  write_csv(rr,file.path(dest,paste0(org,"_",p,".csv")))
  results[[length(results)+1L]] <- rr
  rm(x,fit);gc()
}
write_csv(bind_rows(results),file.path(dest,"registration_sensitivity_results.csv"))
writeLines(c("Separate spells: overlapping or same-day touching registrations combined; positive gaps begin a new spell.",
"Each spell uses its first listing date, address, center, covariates and recomputed preceding-year exposure; current last-terminal-outcome resolution preserved within spells.",
"Model-based variance and a separate person-clustered variance sensitivity are both reported.",
"Concurrency changes at observed registration start/end boundaries. Baseline exposure is unchanged; primary covariates plus concurrent listing count and a zero-open-listing indicator are fitted.",
"The zero-open-listing indicator distinguishes gaps retained by the primary person-organ construction. Maximum future concurrency is never a baseline covariate."),file.path(dest,"methods.txt"))
