source("code/revision_2026/common.R")
base <- baseline_data()
cohort_start <- getOption("revision.recent_start","2020-01-01")
cohort_end <- getOption("revision.recent_end","2022-12-31")
cohort_tag <- paste0(substr(cohort_start,1,4),"_",substr(cohort_end,1,4))
dest <- file.path(revision_dir,getOption("revision.recent_dir","recent_finegray"))
dir.create(dest,showWarnings=FALSE)
results <- list()
for(org in names(organ_names)) for(p in names(pollutants)) {
  term <- pollutants[[p]]
  d <- eligible(base[index_date>=as.Date(cohort_start) & index_date<=as.Date(cohort_end)],org,term)
  d$adverse_event <- as.integer(d$adverse_event==1 & d$followup_days<=3*365.25)
  d$followup_days <- pmin(d$followup_days,3*365.25)
  log_revision("Matched cause-specific ",org," ",p," candidates=",nrow(d))
  warnings <- character()
  fit <- withCallingHandlers(coxph(as.formula(paste("Surv(followup_days,adverse_event) ~",paste(c(term,baseline_adjustment),collapse=" + "))),data=d,ties="efron"),warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w))
  })
  rr <- effect_row(fit,term,paste0("recent_",cohort_tag,"_3yr_cause_specific"),org,p,d)
  rr$convergence_note <- paste(unique(warnings),collapse=" | ")
  results[[length(results)+1L]] <- rr
  write_csv(broom::tidy(fit),file.path(dest,paste0(org,"_",p,"_cause_specific_coefficients.csv")))
}
write_csv(bind_rows(results),file.path(dest,"matched_recent_cause_specific.csv"))
