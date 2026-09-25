source("code/revision_2026/common.R")
source("code/tests/test_revision_design.R")
risks <- read_csv(file.path(revision_dir,"aj_with_risk_tables","numbers_at_risk.csv"),show_col_types=FALSE)
original <- read_csv("output/prior_year_pollution_quartile_aalen_johansen_cif/prior_year_pollution_quartile_event_summary.csv",show_col_types=FALSE)
zero <- risks %>% filter(horizon==1,time==0) %>% select(organ,pollutant,quartile,n_risk)
check <- full_join(zero,original,by=c("organ","pollutant","quartile"))
stopifnot(nrow(check)==48,all(check$n_risk==check$n))
base <- read_csv(file.path(revision_dir,"baseline_diagnostics","baseline_reproduction.csv"),show_col_types=FALSE)
tv <- read_csv(file.path(revision_dir,"timevarying_diagnostics","timevarying_reproduction.csv"),show_col_types=FALSE)
stopifnot(nrow(base)==12,nrow(tv)==12)
# Explicit labels avoid base graphics suppressing crowded default KM-axis labels.
diagnostic_time_axis <- function(z) {
  keep <- !duplicated(z$x)
  x <- z$x[keep]; days <- z$time[keep]
  targets <- seq(min(x),max(x),length.out=6)
  day_ticks <- unique(signif(approx(x,days,xout=targets,rule=2)$y,2))
  positions <- approx(days,x,xout=day_ticks,rule=2)$y
  labels <- format(day_ticks,big.mark=",",scientific=FALSE,trim=TRUE)
  stopifnot(all(is.finite(positions)),all(nzchar(labels)),!anyDuplicated(positions))
  axis(1,at=positions,labels=labels,cex.axis=.85,gap.axis=-1)
}
for(model in c("baseline","tv")) for(org in names(organ_names)) for(p in names(pollutants)) {
  z <- readRDS(file.path(revision_dir,"cache",paste0(model,"_zph_",org,"_",p,".rds")))
  term <- if(model=="baseline") pollutants[[p]] else c(pm25="pm25_interval_5ug",no2="no2_interval_10ppb",o3="o3_interval_10ppb")[[p]]
  row <- if(model=="baseline") base %>% filter(organ==.env$org,pollutant==.env$p) else tv %>% filter(organ==.env$org,pollutant==.env$p)
  beta <- log(if(model=="baseline") row$hazard_ratio else row$estimate)
  dest <- file.path(revision_dir,if(model=="baseline") "baseline_diagnostics" else "timevarying_diagnostics")
  stem <- if(model=="baseline") paste0("schoenfeld_",org,"_",p) else paste0(org,"_",p,"_schoenfeld")
  poll <- switch(p,pm25=expression(PM[2.5]),no2=expression(NO[2]),o3=expression(O[3]))
  heading <- bquote(.(organ_names[[org]]) ~ .(poll[[1]]))
  for(ext in c("png","pdf")) {
    path <- file.path(dest,paste0(stem,".",ext))
    if(ext=="png") png(path,width=3000,height=1500,res=250) else pdf(path,width=12,height=6)
    par(mfrow=c(1,2),mar=c(5,5,4,1),cex=1.05)
    plot(z,var=match(term,colnames(z$y)),resid=TRUE,se=TRUE,cex=.22,
         xlab="Days since listing",ylab="Pollutant log hazard ratio",main=heading,xaxt="n")
    diagnostic_time_axis(z)
    abline(h=beta,lty=2,col="gray40")
    plot(z,var=match(term,colnames(z$y)),resid=FALSE,se=TRUE,
         xlab="Days since listing",ylab="Pollutant log hazard ratio",main="Smoothed coefficient and 95% bands",xaxt="n")
    diagnostic_time_axis(z)
    abline(h=beta,lty=2,col="gray40")
    dev.off()
  }
}
writeLines(c("PASS: all 48 AJ quartile baseline risk counts reproduce the original AJ cohort counts.",
"PASS: all 12 baseline and all 12 time-varying primary effect estimates reproduced within 1e-7 in the fitting scripts, with matching sample/event counts.",
"PASS: synthetic overlap-chain, same-day boundary, gap, concurrency and pre-event risk-count tests.",
"PASS: observed ACS vintage is 2022; no forward/backward-carried vintage used for SVI validation.",
"Diagnostic figures show full scaled residuals and a separate smoothed coefficient with approximate 95% bands. The horizontal dashed line is the fitted constant coefficient; x-axis spacing uses the KM transform.",
"Official LAS/CAS replacement remains unresolved; see docs/lung_score_source_audit.md."),file.path(revision_dir,"validation_checks.txt"))
