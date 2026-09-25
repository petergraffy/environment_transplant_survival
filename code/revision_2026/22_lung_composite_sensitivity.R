source('code/revision_2026/common.R')
source('code/revision_2026/lung_severity_composite.R')
setDTthreads(2)
dest<-file.path(revision_dir,'lung_composite_sensitivity');dir.create(dest,showWarnings=FALSE)
cache<-file.path(revision_dir,'cache','lung_composite');dir.create(cache,showWarnings=FALSE)
transition<-as.Date('2023-03-09')
env<-load_script_prefix('code/79_timevarying_pollution_cox_svi.R','analysis_dat')
allbase<-baseline_data();b<-copy(allbase[WL_ORG=='LU'])
multi<-allbase[,.(multi_organ_candidate=as.integer(uniqueN(WL_ORG)>1)),by=PERS_ID]
b<-multi[b,on='PERS_ID']
b[,`:=`(candidate_zip=env$clean_zip(candidate_zip),listing_era=factor(ifelse(index_date<transition,'Pre-CAS','Post-CAS'),levels=c('Pre-CAS','Post-CAS')),
 age_group=cut(age,c(-Inf,17,39,59,Inf),labels=c('<18','18-39','40-59','60+')),
 race_group=fcase(as.character(race)=='WHITE','White',as.character(race)=='BLACK','Black',as.character(race)=='ASIAN','Asian',default='Other/Unknown'))]
raw<-readRDS(file.path(revision_dir,'cache','lung_las_raw.rds'))
raw<-merge(raw,b[,.(PERS_ID,index_date,observed_end_date,baseline_px=PX_ID)],by='PERS_ID')
raw[,score_date:=as.Date(CAN_LISTING_DT)]
raw<-raw[!is.na(score_date)&score_date>=index_date&score_date<=observed_end_date]
for(v in c('CAN_BMI','CAN_PULM_ART_SYST','CAN_PULM_ART_MEAN'))raw[!is.finite(get(v))|get(v)<=0,(v):=NA_real_]
raw[!CAN_FUNCTN_STAT %in% c(1:3,seq(2010,2100,10),seq(4010,4100,10)),CAN_FUNCTN_STAT:=NA_real_]
raw[,ventilation:=as.numeric(CAN_VENTILATOR)]
raw[!ventilation %in% c(0,1),ventilation:=NA_real_]
raw[is.na(ventilation)&CAN_ON_VENTILATOR %in% c('Y','N'),ventilation:=as.numeric(CAN_ON_VENTILATOR=='Y')]
base_raw<-merge(b[,.(PERS_ID,PX_ID,index_date,listing_era)],raw[,!c('index_date'),with=FALSE],by=c('PERS_ID','PX_ID'),all.x=TRUE)
stopifnot(nrow(base_raw)==nrow(b))
bc<-lung_composite_components(base_raw)
reference<-vapply(bc,function(x)median(x[base_raw$index_date<transition],na.rm=TRUE),numeric(1))
fwrite(data.table(component=names(reference),fixed_pretransition_median=unname(reference)),file.path(dest,'imputation_reference.csv'))
baseline_score<-impute_lung_components(bc,reference)
missing_terms<-grep('^missing_',names(baseline_score),value=TRUE)
base_summary<-rbindlist(lapply(c('All','Pre-CAS','Post-CAS'),function(e){ix<-if(e=='All')rep(TRUE,nrow(bc)) else base_raw$listing_era==e;rbindlist(lapply(names(bc),function(v){x<-bc[[v]][ix];data.table(era=e,component=v,N=length(x),missing_n=sum(!is.finite(x)),missing_percent=100*mean(!is.finite(x)),median=median(x,na.rm=TRUE),q25=quantile(x,.25,na.rm=TRUE),q75=quantile(x,.75,na.rm=TRUE))}))}))
fwrite(base_summary,file.path(dest,'table1_components_missingness.csv'))
score_summary<-rbindlist(lapply(c('All','Pre-CAS','Post-CAS'),function(e){ix<-if(e=='All')rep(TRUE,nrow(bc)) else base_raw$listing_era==e;x<-baseline_score$lung_composite[ix];data.table(era=e,N=length(x),median=median(x),q25=quantile(x,.25),q75=quantile(x,.75))}))
fwrite(score_summary,file.path(dest,'table1_composite_summary.csv'))
raw[,observed_inputs:=rowSums(!is.na(.SD)),.SDcols=c('CAN_DGN','CAN_BMI','CAN_FUNCTN_STAT','ventilation','CAN_PULM_ART_SYST','CAN_PULM_ART_MEAN')]
raw[,baseline_record:=as.integer(PX_ID==baseline_px)]
setorder(raw,PERS_ID,score_date,-baseline_record,-observed_inputs,PX_ID)
raw<-unique(raw,by=c('PERS_ID','score_date'))
carry<-c('CAN_DGN','CAN_BMI','CAN_FUNCTN_STAT','ventilation','CAN_PULM_ART_SYST','CAN_PULM_ART_MEAN')
raw[,(carry):=lapply(.SD,function(x)nafill(as.numeric(x),type='locf')),by=PERS_ID,.SDcols=carry]
scores<-impute_lung_components(lung_composite_components(raw),reference)
events<-cbind(raw[,.(PERS_ID,score_date,index_date,effective_date=score_date+as.integer(score_date>index_date))],as.data.table(scores))
setorder(events,PERS_ID,effective_date)
events[,changed:=c(TRUE,diff(lung_composite)!=0)|rowSums(abs(as.matrix(.SD)-rbind(as.matrix(.SD)[1,,drop=FALSE],head(as.matrix(.SD),-1))))>0,by=PERS_ID,.SDcols=missing_terms]
events<-events[changed==TRUE]
saveRDS(events,file.path(cache,'score_events.rds'))
log_revision('Constructing monthly intervals with composite updates')
d<-env$make_month_intervals(b,as.Date('2025-12-31'))
updates<-unique(rbind(events[,.(PERS_ID,boundary=as.numeric(effective_date-index_date))],b[,.(PERS_ID,boundary=as.numeric(transition-index_date))]))
splits<-updates[d,on=.(PERS_ID,boundary>tstart,boundary<tstop),nomatch=0,allow.cartesian=TRUE,.(interval_row_id=i.interval_row_id,point=x.boundary)]
starts<-unique(rbind(d[,.(interval_row_id,point=tstart)],splits));setorder(starts,interval_row_id,point)
starts<-merge(starts,d[,.(interval_row_id,old_stop=tstop)],by='interval_row_id')
starts[,new_stop:=shift(point,type='lead',fill=old_stop[1]),by=interval_row_id]
d<-merge(d,starts,by='interval_row_id',allow.cartesian=TRUE)
d[,`:=`(tv_adverse_event=as.integer(tv_adverse_event==1&tstop==new_stop),tstart=point,tstop=new_stop)]
d[,`:=`(start_date=index_date+tstart,age_interval=age+tstart/365.25)]
d[,calendar_era:=factor(ifelse(start_date<transition,'Pre-CAS','Post-CAS'),levels=c('Pre-CAS','Post-CAS'))]
d<-merge(d,b[,.(PERS_ID,listing_era)],by='PERS_ID')
setkey(events,PERS_ID,effective_date)
for(v in c('lung_composite',missing_terms))d[,(v):=events[.SD,on=.(PERS_ID,effective_date=start_date),roll=Inf,get(v)]]
d[,score_source_date:=events[.SD,on=.(PERS_ID,effective_date=start_date),roll=Inf,x.score_date]]
stopifnot(!anyNA(d$lung_composite),!any(d$score_source_date>d$start_date))
# Reuse only unchanged monthly pollution values, never LAS/CAS values or offer dates.
previous<-readRDS(file.path(revision_dir,'cache','lung_las_cas','scored_monthly_intervals.rds'))
terms<-c(pm25='pm25_interval_5ug',no2='no2_interval_10ppb',o3='o3_interval_10ppb')
ex<-unique(previous[,c('PERS_ID','year','month',unname(terms)),with=FALSE]);stopifnot(!anyDuplicated(ex[,.(PERS_ID,year,month)]))
d<-merge(d,ex,by=c('PERS_ID','year','month'),all.x=TRUE)
setorder(d,PERS_ID,tstart)
check<-d[,.(start=min(tstart),days=sum(tstop-tstart),end=max(tstop),events=sum(tv_adverse_event),contiguous=all(tstart[-1]==head(tstop,-1))),by=PERS_ID]
stopifnot(all(check$start==0),all(check$days==check$end),all(check$contiguous),all(check$events<=1),sum(d$tv_adverse_event)==sum(previous$tv_adverse_event),sum(d$tstop-d$tstart)==sum(previous$tstop-previous$tstart))
saveRDS(d,file.path(cache,'intervals.rds'));rm(previous);gc()
fwrite(data.table(metric=c('Candidates','Candidates with changed postbaseline score or missingness','Score records','Future score observations used'),value=c(nrow(b),uniqueN(events[score_date>index_date,PERS_ID]),nrow(events),sum(d$score_source_date>d$start_date))),file.path(dest,'score_update_audit.csv'))

results<-list();fullcoef<-list();model_logs<-list();interaction_results<-list()
fit_one<-function(x,exposures,model,subgroup='All',adjust=c('age','sex','race','zcta_svi_proxy'),center=TRUE) {
  log_revision('Fit ',model,' / ',subgroup,' / ',paste(exposures,collapse='+'))
  nuisance<-c(adjust,'lung_composite',missing_terms)
  needed<-unique(c(exposures,nuisance,'listing_center'))
  x<-x[complete.cases(x[,..needed])]
  nuisance<-nuisance[vapply(nuisance,function(v)uniqueN(x[[v]])>1,logical(1))]
  rhs<-c(exposures,nuisance,if(center)'strata(listing_center)')
  f<-as.formula(paste('Surv(tstart,tstop,tv_adverse_event)~',paste(rhs,collapse='+')))
  warnings<-character()
  fit<-tryCatch(withCallingHandlers(coxph(f,data=droplevels(as.data.frame(x)),ties='efron',robust=FALSE),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
  status<-if(inherits(fit,'error'))conditionMessage(fit) else if(length(warnings))paste(unique(warnings),collapse='; ') else 'OK'
  model_logs[[length(model_logs)+1L]]<<-data.table(model,subgroup,exposures=paste(exposures,collapse='+'),status,formula=paste(deparse(f),collapse=' '))
  if(!inherits(fit,'error')) {
    z<-as.data.table(broom::tidy(fit,conf.int=TRUE));z[,`:=`(model=model,subgroup=subgroup,N=uniqueN(x$PERS_ID),events=sum(x$tv_adverse_event),HR=exp(estimate),HR_low=exp(conf.low),HR_high=exp(conf.high),status=status)]
    fullcoef[[length(fullcoef)+1L]]<<-z;results[[length(results)+1L]]<<-z[term %in% exposures]
  }
  fwrite(rbindlist(results),file.path(dest,'pollutant_results.csv'));fwrite(rbindlist(fullcoef),file.path(dest,'full_coefficients.csv'));fwrite(rbindlist(model_logs),file.path(dest,'model_log.csv'))
}
for(p in names(terms)) {
  end<-if(p=='no2')as.Date('2025-12-31') else as.Date('2024-12-31');x<-d[start_date<=end]
  fit_one(x,terms[[p]],'Single pollutant')
  fit_one(x,terms[[p]],'No center strata',center=FALSE)
  fit_one(x[multi_organ_candidate==0],terms[[p]],'Exclude multi-organ')
  fit_one(x,terms[[p]],'Adjust multi-organ',adjust=c('age','sex','race','zcta_svi_proxy','multi_organ_candidate'))
  for(g in levels(b$age_group))fit_one(x[age_group==g],terms[[p]],'Age subgroup',g,c('sex','race_group','zcta_svi_proxy'))
  for(g in c('F','M'))fit_one(x[as.character(sex)==g],terms[[p]],'Sex subgroup',if(g=='F')'Female' else 'Male',c('age_interval','race_group','zcta_svi_proxy'))
  for(g in c('White','Black','Asian'))fit_one(x[race_group==g],terms[[p]],'Race subgroup',g,c('age_interval','sex','zcta_svi_proxy'))
  for(era_var in c('listing_era','calendar_era')) {
    for(e in c('Pre-CAS','Post-CAS'))fit_one(x[get(era_var)==e],terms[[p]],era_var,e)
    need<-c(terms[[p]],'age','sex','race','zcta_svi_proxy','lung_composite',missing_terms,'listing_center',era_var)
    ix<-x[complete.cases(x[,..need])];ix[,era:=get(era_var)]
    adj<-c('age','sex','race','zcta_svi_proxy','lung_composite',missing_terms)
    adj<-adj[vapply(adj,function(v)uniqueN(ix[[v]])>1,logical(1))]
    adj<-c(adj,'strata(listing_center,era)')
    f0<-as.formula(paste('Surv(tstart,tstop,tv_adverse_event)~',paste(c(terms[[p]],adj),collapse='+')))
    f1<-as.formula(paste('Surv(tstart,tstop,tv_adverse_event)~',paste(c(paste0(terms[[p]],':era'),adj),collapse='+')))
    reduced<-coxph(f0,data=as.data.frame(ix),ties='efron',robust=FALSE)
    full<-coxph(f1,data=as.data.frame(ix),ties='efron',robust=FALSE)
    stat<-2*(as.numeric(logLik(full))-as.numeric(logLik(reduced)));df<-attr(logLik(full),'df')-attr(logLik(reduced),'df')
    stopifnot(df==1,stat>=-1e-6)
    interaction_results[[length(interaction_results)+1L]]<-data.table(pollutant=p,era_definition=era_var,chisq=max(stat,0),df,p_value=pchisq(max(stat,0),df,lower.tail=FALSE),N=uniqueN(ix$PERS_ID),events=sum(ix$tv_adverse_event))
    fwrite(rbindlist(interaction_results),file.path(dest,'era_interaction_tests.csv'))
  }
}
fit_one(d[start_date<=as.Date('2024-12-31')],terms[c('pm25','no2')],'PM2.5 + NO2')
fit_one(d[start_date<=as.Date('2024-12-31')],terms,'PM2.5 + NO2 + O3')
log_revision('Composite lung sensitivity models complete')
