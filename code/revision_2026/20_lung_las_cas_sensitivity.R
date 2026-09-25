source('code/revision_2026/common.R')
source('code/revision_2026/las_2021.R')
source('code/saf_paths.R')
setDTthreads(2)
dest<-file.path(revision_dir,'lung_las_cas_sensitivity')
dir.create(dest,showWarnings=FALSE)
cache<-file.path(revision_dir,'cache','lung_las_cas')
dir.create(cache,showWarnings=FALSE)
transition<-as.Date('2023-03-09')
env<-load_script_prefix('code/79_timevarying_pollution_cox_svi.R','analysis_dat')
allbase<-baseline_data()
multi<-allbase[,.(multi_organ_candidate=as.integer(uniqueN(WL_ORG)>1)),by=PERS_ID]
b<-copy(allbase[WL_ORG=='LU'])
b<-multi[b,on='PERS_ID']
b[,candidate_zip:=env$clean_zip(candidate_zip)]
b[,race_group:=fcase(as.character(race)=='WHITE','White',as.character(race)=='BLACK','Black',as.character(race)=='ASIAN','Asian',default='Other/Unknown')]
b[,age_group:=cut(age,c(-Inf,17,39,59,Inf),labels=c('<18','18-39','40-59','60+'))]
stopifnot(!anyDuplicated(b$PERS_ID),!anyDuplicated(b$waitlist_row_id))
raw<-readRDS(file.path(revision_dir,'cache','lung_las_raw.rds'))
raw<-raw[PERS_ID %in% b$PERS_ID]
raw[,score_date:=as.Date(CAN_LISTING_DT)]
raw<-merge(raw,b[,.(PERS_ID,index_date,observed_end_date,baseline_px=PX_ID,baseline_age=age,waitlist_row_id)],by='PERS_ID')
raw<-raw[!is.na(score_date)&score_date>=index_date&score_date<=observed_end_date]
raw[,score_age:=as.numeric(CAN_AGE_AT_LISTING)]
raw[is.na(score_age),score_age:=baseline_age+as.numeric(score_date-index_date)/365.25]
fields<-c('CAN_BMI','CAN_HGT_CM','CAN_WGT_KG','CAN_AT_REST_O2','CAN_SIX_MIN_WALK','CAN_PCO2','CAN_PULM_ART_SYST','CAN_PULM_ART_MEAN','CAN_CARDIAC_OUTPUT')
for(v in fields) raw[!is.finite(get(v))|get(v)<0,(v):=NA_real_]
for(v in c('CAN_BMI','CAN_HGT_CM','CAN_WGT_KG','CAN_CARDIAC_OUTPUT')) raw[get(v)==0,(v):=NA_real_]
raw[CAN_FUNCTN_STAT %in% c(996,998),CAN_FUNCTN_STAT:=NA_real_]
raw[,ventilation:=as.numeric(CAN_VENTILATOR)]
raw[!ventilation %in% c(0,1),ventilation:=NA_real_]
raw[is.na(ventilation)&CAN_ON_VENTILATOR %in% c('Y','N'),ventilation:=as.numeric(CAN_ON_VENTILATOR=='Y')]
raw[,cardiac_index:=CAN_CARDIAC_OUTPUT/sqrt(CAN_HGT_CM*CAN_WGT_KG/3600)]
raw[!is.finite(cardiac_index),cardiac_index:=NA_real_]
raw[,`:=`(dated_creatinine=NA_real_,dated_bilirubin=NA_real_)]
# Baseline Table 1 uses only the first registration, not later measurements.
base_raw<-merge(b[,.(PERS_ID,PX_ID)],raw,by=c('PERS_ID','PX_ID'),all.x=TRUE)
stopifnot(nrow(base_raw)==nrow(b))
base_score<-calculate_las_2021(base_raw)
base_raw[,las:=base_score$las]
base_raw[,diagnosis_group:=base_score$diagnosis_group]
fwrite(raw[, .N,by=.(CAN_DGN,diagnosis_group=las_diagnosis_group(CAN_DGN,CAN_PULM_ART_MEAN))][order(CAN_DGN)],file.path(dest,'diagnosis_mapping_audit.csv'))
# For same-day records prefer the baseline registration, then the most complete form.
raw[,observed_inputs:=rowSums(!is.na(.SD)),.SDcols=c(fields,'CAN_FUNCTN_STAT','ventilation')]
raw[,baseline_record:=as.integer(PX_ID==baseline_px)]
setorder(raw,PERS_ID,score_date,-baseline_record,-observed_inputs,PX_ID)
raw<-unique(raw,by=c('PERS_ID','score_date'))
raw[,new_clinical_information:=observed_inputs>0]
# Carry individual inputs forward only when a subsequent listing supplies new information.
carry<-c(fields,'CAN_DGN','CAN_FUNCTN_STAT','ventilation')
raw[,(carry):=lapply(.SD,function(x) nafill(as.numeric(x),type='locf')),by=PERS_ID,.SDcols=carry]
raw[,cardiac_index:=CAN_CARDIAC_OUTPUT/sqrt(CAN_HGT_CM*CAN_WGT_KG/3600)]
scores<-calculate_las_2021(raw)
raw[,las:=scores$las]
las_events<-raw[(score_date==index_date|new_clinical_information)&score_date<transition,
  .(PERS_ID,score_date,effective_date=score_date+as.integer(score_date>index_date),las)]
las_events<-las_events[is.finite(las)]
saveRDS(las_events,file.path(cache,'las_events.rds'))

casfile<-file.path(cache,'cas_events.rds')
if(file.exists(casfile)) cas_events<-readRDS(casfile) else {
  paths<-get_saf_paths(release='q1_2026')
  files<-list.files(file.path(paths$saf_dir,'ptr_lu_2018_2025'),pattern='202[345]0101.*sas7bdat$',full.names=TRUE)
  stopifnot(length(files)==3)
  links<-unique(raw[,.(PX_ID,PERS_ID)])
  # Include registrations omitted by same-day prioritization in the score history.
  links<-unique(readRDS(file.path(revision_dir,'cache','lung_las_raw.rds'))[PERS_ID %in% b$PERS_ID,.(PX_ID,PERS_ID)])
  stopifnot(!anyDuplicated(links$PX_ID))
  cas_parts<-list()
  for(f in files) {
    log_revision('Reading offer CAS: ',basename(f))
    z<-as.data.table(haven::read_sas(f,col_select=c(PX_ID,MATCH_SUBMIT_DT,PTR_TOT_SCORE)))
    z<-merge(z,links,by='PX_ID')
    z[,score_date:=as.Date(MATCH_SUBMIT_DT)]
    z<-z[score_date>=transition&is.finite(PTR_TOT_SCORE)&PTR_TOT_SCORE>=0&PTR_TOT_SCORE<=100]
    cas_parts[[length(cas_parts)+1L]]<-unique(z[,.(PERS_ID,score_date,PTR_TOT_SCORE)])
    rm(z);gc()
  }
  cas_events<-unique(rbindlist(cas_parts))[,.(cas=median(PTR_TOT_SCORE),distinct_offer_scores=.N),by=.(PERS_ID,score_date)]
  stopifnot(!anyDuplicated(cas_events[,.(PERS_ID,score_date)]))
  cas_events[,effective_date:=score_date+1L]
  cas_events<-merge(cas_events,b[,.(PERS_ID,index_date,observed_end_date)],by='PERS_ID')
  cas_events<-cas_events[effective_date>=index_date&effective_date<=observed_end_date]
  saveRDS(cas_events,casfile)
}

summary_value<-function(x) {
  x<-x[is.finite(x)];if(!length(x))return('Not available')
  sprintf('%.2f (%.2f-%.2f)',median(x),quantile(x,.25),quantile(x,.75))
}
table_fields<-c(age='score_age',diagnosis='CAN_DGN',BMI='CAN_BMI',height_cm='CAN_HGT_CM',weight_kg='CAN_WGT_KG',oxygen_L_min='CAN_AT_REST_O2',walk_feet='CAN_SIX_MIN_WALK',PCO2_mmHg='CAN_PCO2',PA_systolic_mmHg='CAN_PULM_ART_SYST',PA_mean_mmHg='CAN_PULM_ART_MEAN',cardiac_output_L_min='CAN_CARDIAC_OUTPUT',derived_cardiac_index='cardiac_index',functional_status='CAN_FUNCTN_STAT',mechanical_ventilation='ventilation',creatinine_mg_dL='dated_creatinine',bilirubin_mg_dL='dated_bilirubin',reconstructed_LAS='las')
handling<-c('Age at clinical listing form','Explicit diagnosis mapping; unmapped => LAS unavailable','Missing: 100 kg/m2','Used only for derived cardiac index','Used only for derived cardiac index','Missing: WL 0; TX 26.33 L/min','Missing: WL 4000; TX 0 feet','Minimum/default 40 mmHg; no serial rise inferred','Minimum/default 20 mmHg','Missing sarcoidosis pressure: group A','Used with height and weight to derive cardiac index','CO / sqrt(height_cm * weight_kg / 3600); missing: 3','No assistance if code 1 or Karnofsky/Lansky >=70; missing: no assistance','CAN_VENTILATOR; fallback CAN_ON_VENTILATOR; missing: WL 0 / TX 1','No dated LU value: WL 0.1 / TX 40 for adults, no contribution under18','Unavailable: 0.7 mg/dL','2021 formula; unavailable diagnosis or age<12 => indicator')
t1<-rbindlist(lapply(seq_along(table_fields),function(i){x<-base_raw[[table_fields[i]]];data.table(variable=names(table_fields)[i],SRTR_field=table_fields[i],N=nrow(b),observed=sum(is.finite(x)),missing_percent=100*mean(!is.finite(x)),median_IQR=summary_value(x),handling=handling[i])}))
cas_baseline_n<-sum(b$index_date>=transition)
t1<-rbind(t1,data.table(variable='CAS at listing (post-transition candidates)',SRTR_field='PTR_TOT_SCORE',N=cas_baseline_n,observed=0,missing_percent=100,median_IQR='Not available before first eligible offer',handling='No backward filling from later offer; CAS availability indicator until first score'))
fwrite(t1,file.path(dest,'table1_lung_las_cas_inputs.csv'))

interval_file<-file.path(cache,'scored_monthly_intervals.rds')
if(file.exists(interval_file)) intervals<-readRDS(interval_file) else {
  log_revision('Building monthly pollution intervals and score-update boundaries')
  intervals<-env$make_month_intervals(b,as.Date('2025-12-31'))
  le<-merge(las_events,b[,.(PERS_ID,index_date)],by='PERS_ID')
  ce<-copy(cas_events)
  updates<-rbindlist(list(le[,.(PERS_ID,boundary=as.numeric(effective_date-index_date))],ce[,.(PERS_ID,boundary=as.numeric(effective_date-index_date))],b[,.(PERS_ID,boundary=as.numeric(transition-index_date))]))
  updates<-unique(updates)
  # Split each monthly row only at score dates strictly inside that row.
  splits<-updates[intervals,on=.(PERS_ID,boundary>tstart,boundary<tstop),nomatch=0,allow.cartesian=TRUE,
    .(interval_row_id=i.interval_row_id,point=x.boundary)]
  starts<-rbind(intervals[,.(interval_row_id,point=tstart)],splits)
  setorder(starts,interval_row_id,point);starts<-unique(starts)
  starts<-merge(starts,intervals[,.(interval_row_id,original_stop=tstop)],by='interval_row_id')
  starts[,new_stop:=shift(point,type='lead',fill=original_stop[1]),by=interval_row_id]
  intervals<-merge(intervals,starts,by='interval_row_id',allow.cartesian=TRUE)
  intervals[,`:=`(tstart=point,tstop=new_stop,tv_adverse_event=as.integer(tv_adverse_event==1&tstop==new_stop))]
  intervals[,start_date:=index_date+tstart]
  setkey(las_events,PERS_ID,effective_date);setkey(cas_events,PERS_ID,effective_date)
  intervals[,las_tv:=las_events[.SD,on=.(PERS_ID,effective_date=start_date),roll=Inf,x.las]]
  intervals[,cas_tv:=cas_events[.SD,on=.(PERS_ID,effective_date=start_date),roll=Inf,x.cas]]
  intervals[,post_cas:=as.integer(start_date>=transition)]
  intervals[,`:=`(las_unavailable=as.integer(post_cas==0&is.na(las_tv)),cas_unavailable=as.integer(post_cas==1&is.na(cas_tv)),
    las_10=ifelse(post_cas==0&is.finite(las_tv),las_tv/10,0),cas_10=ifelse(post_cas==1&is.finite(cas_tv),cas_tv/10,0),age_interval=age+tstart/365.25)]
  ex<-env$read_exposure_tables()
  for(nm in c('pm25','o3','no2_monthly')) {
    value<-c(pm25='pm25_interval_ug_m3',o3='o3_interval_ppb',no2_monthly='no2_interval_ppb')[[nm]]
    intervals[,(value):=ex[[nm]][.SD,on=.(zip=candidate_zip,year,month),get(value)]]
  }
  intervals[,no2_annual_ppb:=ex$no2_annual[.SD,on=.(zip=candidate_zip,year),no2_annual_ppb]]
  intervals[,`:=`(pm25_interval_5ug=pm25_interval_ug_m3/5,o3_interval_10ppb=o3_interval_ppb/10,no2_interval_10ppb=fcoalesce(no2_interval_ppb,no2_annual_ppb)/10)]
  saveRDS(intervals,interval_file)
}
setorder(intervals,PERS_ID,tstart)
stopifnot(all(intervals$tstop>intervals$tstart))
checks<-intervals[,.(start=min(tstart),days=sum(tstop-tstart),end=max(tstop),events=sum(tv_adverse_event),contiguous=all(tstart[-1]==head(tstop,-1))),by=PERS_ID]
stopifnot(all(checks$contiguous),all(checks$start==0),all(checks$events<=1),all(checks$days==checks$end))
original<-env$make_month_intervals(b,as.Date('2025-12-31'))
stopifnot(sum(original$tv_adverse_event)==sum(intervals$tv_adverse_event),sum(original$tstop-original$tstart)==sum(intervals$tstop-intervals$tstart))
coverage<-intervals[,.(candidates=uniqueN(PERS_ID),intervals=.N,person_days=sum(tstop-tstart),LAS_available_days=sum((tstop-tstart)*(post_cas==0&las_unavailable==0)),CAS_available_days=sum((tstop-tstart)*(post_cas==1&cas_unavailable==0)),events=sum(tv_adverse_event)),by=post_cas]
fwrite(coverage,file.path(dest,'score_coverage.csv'))
updates_summary<-data.table(metric=c('lung_candidates','candidates_with_reconstructed_baseline_LAS','candidates_with_postbaseline_LAS_record','candidates_with_eligible_CAS_record','LAS_records','CAS_candidate_days'),value=c(nrow(b),sum(is.finite(base_raw$las)),uniqueN(raw[score_date>index_date&score_date<transition&is.finite(las),PERS_ID]),uniqueN(cas_events$PERS_ID),nrow(las_events),nrow(cas_events)))
fwrite(updates_summary,file.path(dest,'score_updates.csv'))

results<-list();coefficients<-list();logs<-list()
fit_one<-function(d,terms,model,subgroup='All',adjust=c('age','sex','race','zcta_svi_proxy'),center=TRUE) {
  log_revision('Fit ',model,' / ',subgroup,' / ',paste(terms,collapse='+'))
  severity<-c('post_cas','las_10','cas_10','las_unavailable','cas_unavailable')
  need<-unique(c(terms,adjust,severity,'listing_center'))
  d<-d[complete.cases(d[,..need])]
  # Remove constant nuisance terms in small subgroups, retaining exposure terms.
  nuisance<-c(adjust,severity)
  nuisance<-nuisance[vapply(nuisance,function(v)uniqueN(d[[v]])>1,logical(1))]
  rhs<-c(terms,nuisance,if(center)'strata(listing_center)')
  formula<-as.formula(paste('Surv(tstart,tstop,tv_adverse_event)~',paste(rhs,collapse='+')))
  warn<-character()
  if(sum(d$tv_adverse_event)<50||uniqueN(d$listing_center)<2) {
    logs[[length(logs)+1L]]<<-data.table(model,subgroup,status='Not estimated: fewer than 50 events or 2 centers',formula=paste(deparse(formula),collapse=' '));return(invisible(NULL))
  }
  fit<-tryCatch(withCallingHandlers(coxph(formula,data=droplevels(as.data.frame(d)),ties='efron',robust=FALSE,control=coxph.control(iter.max=50)),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
  if(inherits(fit,'error')) {logs[[length(logs)+1L]]<<-data.table(model,subgroup,status=conditionMessage(fit),formula=paste(deparse(formula),collapse=' '));return(invisible(NULL))}
  z<-as.data.table(broom::tidy(fit,conf.int=TRUE))
  z[,`:=`(model=model,subgroup=subgroup,candidates=uniqueN(d$PERS_ID),events=sum(d$tv_adverse_event),HR=exp(estimate),HR_low=exp(conf.low),HR_high=exp(conf.high))]
  coefficients[[length(coefficients)+1L]]<<-z
  results[[length(results)+1L]]<<-z[term %in% terms]
  logs[[length(logs)+1L]]<<-data.table(model,subgroup,status=if(length(warn))paste(unique(warn),collapse='; ') else 'OK',formula=paste(deparse(formula),collapse=' '))
  fwrite(rbindlist(results),file.path(dest,'pollutant_results.csv'))
  fwrite(rbindlist(coefficients),file.path(dest,'full_coefficients.csv'))
  fwrite(rbindlist(logs),file.path(dest,'model_log.csv'))
  invisible(gc())
}
terms<-c(pm25='pm25_interval_5ug',no2='no2_interval_10ppb',o3='o3_interval_10ppb')
for(p in names(terms)) {
  end<-if(p=='no2')as.Date('2025-12-31') else as.Date('2024-12-31')
  d<-intervals[start_date<=end]
  fit_one(d,terms[[p]],'Single pollutant')
  fit_one(d,terms[[p]],'No center strata',center=FALSE)
  fit_one(d[multi_organ_candidate==0],terms[[p]],'Exclude multi-organ')
  fit_one(d,terms[[p]],'Adjust multi-organ',adjust=c('age','sex','race','zcta_svi_proxy','multi_organ_candidate'))
  for(g in levels(b$age_group))fit_one(d[age_group==g],terms[[p]],'Age subgroup',g,c('sex','race_group','zcta_svi_proxy'))
  for(g in c('F','M'))fit_one(d[as.character(sex)==g],terms[[p]],'Sex subgroup',if(g=='F')'Female' else 'Male',c('age_interval','race_group','zcta_svi_proxy'))
  for(g in c('White','Black','Asian'))fit_one(d[race_group==g],terms[[p]],'Race subgroup',g,c('age_interval','sex','zcta_svi_proxy'))
}
fit_one(intervals[start_date<=as.Date('2024-12-31')],terms[c('pm25','no2')],'PM2.5 + NO2')
fit_one(intervals[start_date<=as.Date('2024-12-31')],terms,'PM2.5 + NO2 + O3')
log_revision('Sensitivity models complete: ',dest)
