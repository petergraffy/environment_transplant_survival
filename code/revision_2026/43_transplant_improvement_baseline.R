source("code/revision_2026/common.R")
setDTthreads(2)
dest <- file.path(revision_dir, "transplant_improvement_baseline")
dir.create(dest, recursive=TRUE, showWarnings=FALSE)
d <- baseline_data()
lung <- readRDS(file.path(revision_dir,"cache","las_extended_followup","prepared.rds"))$cohort
d <- rbindlist(list(d[WL_ORG!="LU"],lung),fill=TRUE)
stopifnot(nrow(d)==1007801L, !anyDuplicated(d[,.(PERS_ID,WL_ORG)]),
  all(d[WL_ORG=="LU",index_date]<=as.Date("2023-03-08")))
components <- readRDS(file.path(revision_dir,"cache","primary_outcome_components.rds"))
d[, component := NULL]
d <- merge(d,components[,.(PERS_ID,WL_ORG,component)],by=c("PERS_ID","WL_ORG"),all.x=TRUE)
stopifnot(!anyNA(d$component),!anyNA(d$transplant_or_improvement),
  all(d$transplant_or_improvement==as.integer(d$component %in% c("transplant","improvement"))),
  all(d$adverse_event==as.integer(d$component %in% c("death","deterioration"))),
  all(d$transplant_or_improvement+d$adverse_event<=1))
fwrite(d[,.(N=.N),by=.(WL_ORG,component)],file.path(dest,"cohort_outcome_components.csv"))
results <- list(); coefficients <- list(); logs <- list()
reference <- fread("output/prior_year_pollution_cox_svi/prior_year_pollution_cox_svi_results.csv")
for(org in names(organ_names)) for(p in names(pollutants)) {
  term <- pollutants[[p]]
  x <- eligible(d,org,term)
  stopifnot(!anyNA(x$transplant_or_improvement))
  form <- as.formula(paste("Surv(followup_days, transplant_or_improvement) ~",
    paste(c(term,baseline_adjustment),collapse=" + ")))
  log_revision(organ_names[[org]]," / ",p," / N=",nrow(x)," / transplant-improvement events=",sum(x$transplant_or_improvement))
  warnings <- character()
  fit <- withCallingHandlers(coxph(form,data=x,ties="efron",robust=FALSE,x=FALSE,y=FALSE),
    warning=function(w) {warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")})
  z <- as.data.table(broom::tidy(fit,conf.int=TRUE))
  z[, `:=`(organ=org,organ_label=organ_names[[org]],pollutant=p,N=nrow(x),
    events=sum(x$transplant_or_improvement),transplants=sum(x$component=="transplant"),
    improvements=sum(x$component=="improvement"),adverse_competing_events=sum(x$adverse_event),
    centers=length(unique(x$listing_center)),formula=paste(deparse(form),collapse=" "),
    status=if(length(warnings)) paste(unique(warnings),collapse="; ") else "OK",
    HR=exp(estimate),HR_low=exp(conf.low),HR_high=exp(conf.high))]
  coefficients[[length(coefficients)+1L]] <- z
  results[[length(results)+1L]] <- z[term==pollutants[[p]]]
  logs[[length(logs)+1L]] <- unique(z[,.(organ,pollutant,N,events,centers,formula,status)])
  fwrite(rbindlist(results),file.path(dest,"pollutant_results.csv"))
  fwrite(rbindlist(coefficients),file.path(dest,"full_coefficients.csv"))
  fwrite(rbindlist(logs),file.path(dest,"model_log.csv"))
  rm(fit,x); gc(verbose=FALSE)
}
r <- rbindlist(results)
stopifnot(nrow(r)==12L,all(is.finite(r$HR)),all(r$events==r$transplants+r$improvements))
lung_reference <- fread(file.path(revision_dir,"las_extended_followup","pollutant_results.csv"))[model=="Baseline" & variant=="Primary"]
lung_reference[, pollutant := names(pollutants)[match(term,pollutants)]]
expected <- rbindlist(list(reference[organ!="LU",.(organ,pollutant=exposure,primary_N=n,primary_adverse_events=adverse_events)],
  lung_reference[,.(organ="LU",pollutant,primary_N=N,primary_adverse_events=events)]))
check <- merge(r[,.(organ,pollutant,N,adverse_competing_events)],expected,by=c("organ","pollutant"))
stopifnot(nrow(check)==12L,all(check$N==check$primary_N),all(check$adverse_competing_events==check$primary_adverse_events))
fwrite(check,file.path(dest,"primary_cohort_reconciliation.csv"))
notes <- c(
  "Separate cause-specific Cox models for heart, kidney, liver and lung, one pollutant per model.",
  "Outcome: transplantation or removal for improvement, using the primary cohort's terminal outcome coding.",
  "Death/deterioration are competing events censored at occurrence; other exits and administrative end are censored as in the primary analysis.",
  "Exposure: preceding 365 days for PM2.5/O3; prior 12 complete calendar months for NO2, with annual approximation before monthly availability.",
  "Adjustment: age, sex, race, ACS vulnerability proxy (zcta_svi_proxy), and strata(listing_center). No listing year or organ score.",
  "PM2.5 HR per 5 ug/m3; NO2 and O3 HR per 10 ppb. Efron ties; model-based standard errors, no robust person clustering.",
  "The same complete-case eligibility as the primary baseline models is used. Nonlung cohorts unchanged; lung entry through March 8, 2023 with extended follow-up up to December 31, 2025.",
  "HR >1 indicates a higher instantaneous rate of transplantation/improvement among candidates still at risk, not a higher absolute probability or a causal benefit.",
  "The primary transplant category is retained, including removal codes classified as transplantation despite operative death; this analysis does not redefine the endpoint as successful post-transplant survival.")
writeLines(notes,file.path(dest,"methods.txt"))
display <- r[,.(Organ=organ_label,Pollutant=pollutant,N,Transplants=transplants,Improvements=improvements,
  `HR (95% CI)`=sprintf("%.2f (%.2f-%.2f)",HR,HR_low,HR_high),
  `P value`=ifelse(p.value<.001,"<.001",sprintf("%.3f",p.value)))]
fwrite(display,file.path(dest,"supplement_table.csv"))
escape <- function(s) gsub(">","&gt;",gsub("<","&lt;",gsub("&","&amp;",as.character(s),fixed=TRUE),fixed=TRUE),fixed=TRUE)
rows <- apply(as.data.frame(display),1,function(v) paste0("<tr>",paste0("<td>",escape(v),"</td>",collapse=""),"</tr>"))
writeLines(c("<!doctype html><meta charset='utf-8'><title>Baseline pollution and transplantation or improvement</title>",
  "<style>body{font:16px Arial;max-width:1150px;margin:35px auto;color:#111}table{border-collapse:collapse;width:100%}th,td{padding:10px;border-bottom:1px solid #bbb;text-align:left}th{background:#eee}</style>",
  "<h1>Baseline air pollution and transplantation or improvement</h1>",
  "<table><thead><tr>",paste0("<th>",escape(names(display)),"</th>",collapse=""),"</tr></thead><tbody>",rows,"</tbody></table>",
  paste0("<p>",escape(notes),"</p>"),"<h2>Model warnings</h2>",
  if(all(r$status=="OK")) "<p>All 12 models completed without warnings.</p>" else paste0("<p>",escape(r[status!="OK",paste(organ,pollutant,status)]),"</p>")),
  file.path(dest,"supplement_table.html"))
print(display); print(r[,.(organ,pollutant,N,events,status)])
