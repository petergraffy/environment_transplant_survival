source('code/revision_2026/common.R')
setDTthreads(2)
cache<-file.path(revision_dir,'cache','lung_las_cas')
d<-readRDS(file.path(cache,'scored_monthly_intervals.rds'))
le<-readRDS(file.path(cache,'las_events.rds'))
ce<-readRDS(file.path(cache,'cas_events.rds'))
setkey(le,PERS_ID,effective_date);setkey(ce,PERS_ID,effective_date)
las_dates<-le[d,on=.(PERS_ID,effective_date=start_date),roll=Inf,x.score_date]
cas_dates<-ce[d,on=.(PERS_ID,effective_date=start_date),roll=Inf,x.score_date]
stopifnot(!any(las_dates>d$start_date,na.rm=TRUE),!any(cas_dates>=d$start_date,na.rm=TRUE))
stopifnot(all(d[post_cas==1,las_10]==0),all(d[post_cas==0,cas_10]==0),
  all(d[las_unavailable==1,las_10]==0),all(d[cas_unavailable==1,cas_10]==0))
z<-fread(file.path(revision_dir,'lung_las_cas_sensitivity','supplement_results.csv'))
stopifnot(nrow(z)==44,!anyDuplicated(z[,.(model,subgroup,term)]))
fwrite(data.table(check=c('Future LAS measurements used','Same-day or future CAS measurements used','Pollutant estimates exported'),value=c(sum(las_dates>d$start_date,na.rm=TRUE),sum(cas_dates>=d$start_date,na.rm=TRUE),nrow(z))),file.path(revision_dir,'lung_las_cas_sensitivity','timing_validation.csv'))
cat('No-lookahead, era separation and complete result-set tests passed.\n')
