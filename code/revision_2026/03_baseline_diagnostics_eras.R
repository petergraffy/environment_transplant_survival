source("code/revision_2026/common.R")
dat <- baseline_data()
dest <- file.path(revision_dir,"baseline_diagnostics")
dir.create(dest,showWarnings=FALSE)
era_cuts <- list(HR=as.Date("2018-10-18"),
  KI=as.Date(c("2014-12-04","2021-03-15")),
  LI=as.Date(c("2016-01-11","2020-02-04","2023-07-13")),
  LU=as.Date(c("2017-11-24","2023-03-09")))
write_csv(bind_rows(lapply(names(era_cuts),function(org) data.frame(organ=org,cut_date=era_cuts[[org]]))),
          file.path(dest,"policy_era_cutpoints.csv"))
results <- list(); tests <- list(); era_results <- list(); interaction_tests <- list(); bands <- list()
reference <- read_csv("output/prior_year_pollution_cox_svi/prior_year_pollution_cox_svi_results.csv",show_col_types=FALSE)
for(org in names(organ_names)) for(p in names(pollutants)) {
  term <- pollutants[[p]]; d <- eligible(dat,org,term)
  log_revision("Baseline PH/era audit ",org," ",p," n=",nrow(d))
  form <- as.formula(paste("Surv(followup_days,adverse_event) ~",paste(c(term,baseline_adjustment),collapse=" + ")))
  fit_context <- "primary_baseline"
  fit <- coxph(form,d,ties="efron",x=TRUE,y=TRUE)
  row <- effect_row(fit,term,"baseline_reference_reproduction",org,p,d)
  old <- reference %>% filter(organ==org,exposure==p)
  stopifnot(nrow(old)==1L,abs(log(row$hazard_ratio/old$hazard_ratio))<1e-7,row$n==old$n)
  results[[length(results)+1L]] <- row
  z <- cox.zph(fit,transform="km",terms=TRUE)
  tests[[length(tests)+1L]] <- data.frame(organ=org,pollutant=p,term=rownames(z$table),z$table,row.names=NULL)
  saveRDS(z,file.path(revision_dir,"cache",paste0("baseline_zph_",org,"_",p,".rds")))
  for(ext in c("png","pdf")) {
    path <- file.path(dest,paste0("schoenfeld_",org,"_",p,".",ext))
    if(ext=="png") png(path,width=1600,height=1100,res=180) else pdf(path,width=9,height=6)
    plot(z,var=match(term,colnames(z$y)),resid=TRUE,se=TRUE,
         xlab="Days since listing",ylab="Pollutant log hazard ratio",main=paste(organ_names[[org]],p))
    abline(h=coef(fit)[term],lty=2,col="grey40"); dev.off()
  }
  # Estimate time-band-specific coefficients without asserting PH from P values alone.
  bounds <- c(0,365.25,3*365.25,5*365.25,Inf)
  for(k in 1:4) {
    dd <- d[d$followup_days>bounds[k],]
    dd$start <- bounds[k]; dd$stop <- pmin(dd$followup_days,bounds[k+1])
    dd$band_event <- as.integer(dd$adverse_event==1 & dd$followup_days<=bounds[k+1])
    if(sum(dd$band_event)<20) next
    fit_context <- paste0("baseline_band_",k)
    f <- coxph(as.formula(paste("Surv(start,stop,band_event) ~",paste(c(term,baseline_adjustment),collapse=" + "))),dd,ties="efron")
    rr <- effect_row(f,term,paste0("baseline_followup_band_",k),org,p,dd)
    rr$events <- sum(dd$band_event); rr$lower_year <- bounds[k]/365.25; rr$upper_year <- bounds[k+1]/365.25
    bands[[length(bands)+1L]] <- rr
  }
  # Adult heart allocation policy does not define pediatric status transitions.
  de <- if(org=="HR") d[d$age>=18,] else d
  cuts <- c(as.Date("1900-01-01"),era_cuts[[org]],as.Date("2100-01-01"))
  de$era <- cut(de$index_date,breaks=cuts,right=FALSE,labels=paste0("Era ",seq_len(length(cuts)-1)))
  for(e in levels(de$era)) {
    dd <- droplevels(de[de$era==e,]); if(sum(dd$adverse_event)<20) next
    fit_context <- paste0("listing_",e)
    f <- coxph(form,dd,ties="efron")
    rr <- effect_row(f,term,paste0("listing_",e),org,p,dd)
    rr$adult_heart_only <- org=="HR"
    era_results[[length(era_results)+1L]] <- rr
  }
  adj <- c("age","sex","race","zcta_svi_proxy","strata(listing_center,era)")
  fit_context <- "era_interaction_reduced"
  reduced <- coxph(as.formula(paste("Surv(followup_days,adverse_event) ~",paste(c(term,adj),collapse=" + "))),de)
  fit_context <- "era_interaction_full"
  full <- coxph(as.formula(paste("Surv(followup_days,adverse_event) ~",paste(c(paste0(term,":era"),adj),collapse=" + "))),de)
  df <- sum(!is.na(coef(full)))-sum(!is.na(coef(reduced)))
  statistic <- 2*(full$loglik[2]-reduced$loglik[2])
  interaction_tests[[length(interaction_tests)+1L]] <- data.frame(organ=org,pollutant=p,chisq=statistic,df=df,p_value=pchisq(statistic,df,lower.tail=FALSE),adult_heart_only=org=="HR")
  rm(fit,z,full,reduced); gc()
  write_csv(bind_rows(results),file.path(dest,"baseline_reproduction.csv"))
  write_csv(bind_rows(tests),file.path(dest,"schoenfeld_tests.csv"))
  write_csv(bind_rows(bands),file.path(dest,"followup_band_effects.csv"))
  write_csv(bind_rows(era_results),file.path(dest,"listing_era_effects.csv"))
  write_csv(bind_rows(interaction_tests),file.path(dest,"pollutant_by_era_tests.csv"))
}
