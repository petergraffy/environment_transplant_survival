source('code/revision_2026/common.R')
setDTthreads(2)
dest<-file.path(revision_dir,'aj_risk_audit');dir.create(dest,showWarnings=FALSE)
base<-baseline_data()
fields<-c('WL_ORG','PERS_ID','index_date','observed_end_date','followup_days','adverse_event','transplant_or_improvement','other_exit','event_type','pm25_prior_ug_m3','no2_prior_ppb','o3_prior_ppb')
base<-copy(base[,..fields])
stopifnot(all(base$followup_days==as.numeric(base$observed_end_date-base$index_date),na.rm=TRUE))
existing<-fread(file.path(revision_dir,'aj_with_risk_tables','numbers_at_risk.csv'))
existing<-unique(existing[,.(organ,pollutant,quartile,time,n_risk)])
saved_cif<-fread('output/prior_year_pollution_quartile_aalen_johansen_cif/prior_year_pollution_quartile_cif_at_1_3_5_10_years.csv')
end_date<-max(base$observed_end_date,na.rm=TRUE)
breakdown<-list();mix<-list();cifs<-list();risk_checks<-list();statuses<-list()
raws<-c('PM2.5'='pm25_prior_ug_m3','NO2'='no2_prior_ppb','O3'='o3_prior_ppb')
qs<-c('Q1 lowest','Q2','Q3','Q4 highest')
for(p in names(raws)) {
  log_revision('Audit ',p)
  col<-raws[[p]]
  d<-copy(base[is.finite(get(col))&followup_days>0&!is.na(adverse_event)&!is.na(transplant_or_improvement)])
  cuts<-quantile(d[[col]],seq(0,1,.25),type=7);cuts[c(1,5)]<-c(-Inf,Inf)
  d[,`:=`(quartile=cut(get(col),cuts,labels=qs,include.lowest=TRUE),organ=unname(organ_names[WL_ORG]),listing_year=as.integer(format(index_date,'%Y')),
    years=followup_days/365.25,status=factor(ifelse(adverse_event==1,'adverse',ifelse(transplant_or_improvement==1,'competing','censor')),levels=c('censor','adverse','competing')))]
  mix[[p]]<-d[,.(N=.N,listing_year_median=as.numeric(median(listing_year)),listing_year_q25=quantile(listing_year,.25),listing_year_q75=quantile(listing_year,.75),listed_2020plus_percent=100*mean(listing_year>=2020),potential_3yr_percent=100*mean(as.numeric(end_date-index_date)>=3*365.25),median_observed_days=as.numeric(median(followup_days))),by=.(organ,quartile)][,pollutant:=p]
  statuses[[p]]<-d[, .N,by=.(organ,quartile,event_type)][,pollutant:=p]
  for(t in c(1,3,5,10)) {
    b<-d[,.(N=.N,at_risk=sum(years>=t),adverse_before=sum(years<t&status=='adverse'),competing_before=sum(years<t&status=='competing'),censored_before=sum(years<t&status=='censor'),
      censored_at_data_end_before=sum(years<t&status=='censor'&observed_end_date==end_date),other_exit_censored_before=sum(years<t&status=='censor'&other_exit==1,na.rm=TRUE)),by=.(organ,quartile)]
    b[,`:=`(pollutant=p,time=t,retained_percent=100*at_risk/N)]
    stopifnot(all(b$N==b$at_risk+b$adverse_before+b$competing_before+b$censored_before))
    breakdown[[length(breakdown)+1L]]<-b
  }
  for(org in unique(d$organ)) {
    x<-d[organ==org]
    fit<-survfit(Surv(years,status)~quartile,data=as.data.frame(x))
    tt<-sort(unique(existing[organ==org&pollutant==p,time]))
    z<-summary(fit,times=tt,extend=TRUE)
    nr<-if(is.matrix(z$n.risk))z$n.risk[,1] else z$n.risk
    audit<-data.table(organ=org,pollutant=p,quartile=sub('^quartile=','',as.character(z$strata)),time=z$time,independent_n_risk=as.numeric(nr),
      adverse_cif=as.numeric(z$pstate[,match('adverse',z$states)]),competing_cif=as.numeric(z$pstate[,match('competing',z$states)]))
    risk_checks[[length(risk_checks)+1L]]<-merge(audit,existing,by=c('organ','pollutant','quartile','time'))
    cifs[[length(cifs)+1L]]<-audit[time %in% c(1,3,5,10)]
  }
}
r<-rbindlist(risk_checks);stopifnot(nrow(r)==nrow(existing),all(r$independent_n_risk==r$n_risk))
c<-rbindlist(cifs)
v<-merge(c,saved_cif,by=c('organ','pollutant','quartile','time'))
stopifnot(nrow(v)==192,all(abs(v$adverse_cif-v$cif_adverse)<1e-10),all(abs(v$competing_cif-v$cif_transplant_or_improvement)<1e-10))
out<-merge(rbindlist(breakdown),c[,.(organ,pollutant,quartile,time,adverse_cif,competing_cif)],by=c('organ','pollutant','quartile','time'))
fwrite(out,file.path(dest,'risk_set_decomposition.csv'))
fwrite(rbindlist(mix),file.path(dest,'listing_year_and_followup_mix.csv'))
fwrite(r,file.path(dest,'independent_risk_count_validation.csv'))
fwrite(rbindlist(statuses),file.path(dest,'event_type_counts.csv'))
summary<-merge(out[time==3&quartile %in% qs[c(1,4)]],rbindlist(mix)[,.(organ,pollutant,quartile,listing_year_median,listed_2020plus_percent,potential_3yr_percent)],by=c('organ','pollutant','quartile'))
fwrite(summary,file.path(dest,'q1_q4_three_year_summary.csv'))
writeLines(c(paste('Latest observed end date:',end_date),'All displayed risk counts reproduced using independent survival::survfit multistate fits.','All saved CIFs at years 1, 3, 5 and 10 reproduced within 1e-10.','At each horizon: starting N = at risk + prior adverse events + prior competing events + prior censoring.','Censoring at the latest observed end date is reported as such, not assumed to prove administrative censoring.'),file.path(dest,'validation.txt'))
print(summary[pollutant=='NO2'&organ %in% c('Heart','Liver','Lung')])
