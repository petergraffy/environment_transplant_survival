source('code/revision_2026/common.R')
source('code/revision_2026/lung_severity_composite.R')
setDTthreads(2)
dest<-file.path(revision_dir,'lung_composite_sensitivity')
d<-readRDS(file.path(revision_dir,'cache','lung_composite','intervals.rds'))
r<-fread(file.path(dest,'pollutant_results.csv'));log<-fread(file.path(dest,'model_log.csv'));tests<-fread(file.path(dest,'era_interaction_tests.csv'))
stopifnot(nrow(log)==53,nrow(tests)==6,!anyDuplicated(r[,.(model,subgroup,term)]))
failed<-log[grepl('overflow',status,fixed=TRUE)]
if(nrow(failed))r<-rbind(r,failed[,.(model,subgroup,term=exposures,status)],fill=TRUE)
stopifnot(nrow(r)==56)
term_map<-c(pm25_interval_5ug='PM2.5',no2_interval_10ppb='NO2',o3_interval_10ppb='O3')
r[,`:=`(Pollutant=unname(term_map[term]),HR_CI=sprintf('%.2f (%.2f-%.2f)',HR,HR_low,HR_high),P=ifelse(p.value<.001,'<.001',sprintf('%.3f',p.value)))]
fwrite(r,file.path(dest,'supplement_results.csv'));fwrite(r,file.path(dest,'supplement_results.tsv'),sep='\t')
comparison<-fread('output/timevarying_pollution_cox_svi/timevarying_pollution_cox_svi_results.csv')[organ=='LU']
validation<-list()
for(p in c('pm25','no2','o3')) {
  term<-c(pm25='pm25_interval_5ug',no2='no2_interval_10ppb',o3='o3_interval_10ppb')[[p]]
  end<-if(p=='no2')as.Date('2025-12-31') else as.Date('2024-12-31')
  x<-d[start_date<=end];needed<-c(term,'age','sex','race','zcta_svi_proxy','organ_score','listing_center')
  x<-x[complete.cases(x[,..needed])]
  fit<-coxph(as.formula(paste('Surv(tstart,tstop,tv_adverse_event)~',paste(c(term,'age','sex','race','zcta_svi_proxy','organ_score','strata(listing_center)'),collapse='+'))),data=as.data.frame(x),ties='efron',robust=FALSE)
  hr<-exp(coef(fit)[term]);old<-comparison[pollutant==p]
  stopifnot(abs(hr-old$hazard_ratio)<1e-6,uniqueN(x$PERS_ID)==old$candidate_organ_episodes,sum(x$tv_adverse_event)==old$adverse_events)
  validation[[p]]<-data.table(pollutant=p,original_HR=old$hazard_ratio,replicated_HR=unname(hr),N=uniqueN(x$PERS_ID),events=sum(x$tv_adverse_event))
}
fwrite(rbindlist(validation),file.path(dest,'original_model_replication.csv'))
era_effects<-list();era_coefficients<-list()
missing_terms<-grep('^missing_',names(d),value=TRUE)
for(p in c('pm25','no2','o3'))for(era_var in c('listing_era','calendar_era')) {
  exposure_term<-c(pm25='pm25_interval_5ug',no2='no2_interval_10ppb',o3='o3_interval_10ppb')[[p]]
  end<-if(p=='no2')as.Date('2025-12-31') else as.Date('2024-12-31')
  x<-copy(d[start_date<=end]);x[,era:=get(era_var)]
  need<-c(exposure_term,'age','sex','race','zcta_svi_proxy','lung_composite',missing_terms,'listing_center','era')
  x<-x[complete.cases(x[,..need])]
  adj<-c('age','sex','race','zcta_svi_proxy','lung_composite',missing_terms)
  adj<-adj[vapply(adj,function(v)uniqueN(x[[v]])>1,logical(1))]
  f<-as.formula(paste('Surv(tstart,tstop,tv_adverse_event)~',paste(c(paste0(exposure_term,':era'),adj,'strata(listing_center,era)'),collapse='+')))
  w<-character()
  fit<-withCallingHandlers(coxph(f,data=as.data.frame(x),ties='efron',robust=FALSE),warning=function(z){w<<-c(w,conditionMessage(z));invokeRestart('muffleWarning')})
  z<-as.data.table(broom::tidy(fit,conf.int=TRUE,exponentiate=TRUE));z[,`:=`(pollutant=p,era_definition=era_var,status=if(length(w))paste(w,collapse='; ') else 'OK')]
  era_coefficients[[length(era_coefficients)+1L]]<-z
  z<-z[startsWith(term,paste0(exposure_term,':era'))]
  stopifnot(nrow(z)==2)
  z[,era:=sub('.*:era','',term)]
  z[,`HR (95% CI)`:=sprintf('%.2f (%.2f-%.2f)',estimate,conf.low,conf.high)]
  z[,interaction_P:=tests[pollutant==p&era_definition==era_var,p_value]]
  era_effects[[length(era_effects)+1L]]<-z
}
era_effects<-rbindlist(era_effects)
fwrite(era_effects,file.path(dest,'era_effects_from_interaction_models.csv'))
fwrite(rbindlist(era_coefficients),file.path(dest,'era_model_full_coefficients.csv'))
coverage<-d[,.(candidates=uniqueN(PERS_ID),person_days=sum(tstop-tstart),events=sum(tv_adverse_event)),by=calendar_era]
fwrite(coverage,file.path(dest,'calendar_era_coverage.csv'))
missing<-fread(file.path(dest,'table1_components_missingness.csv'))
medians<-fread(file.path(dest,'imputation_reference.csv'))
weights<-data.table(component=c('diagnosis','functional_dependence','ventilation','low_BMI','pulmonary_pressure'),
 SRTR_fields=c('CAN_DGN; CAN_PULM_ART_MEAN (sarcoidosis)','CAN_FUNCTN_STAT','CAN_VENTILATOR; CAN_ON_VENTILATOR fallback','CAN_BMI','CAN_PULM_ART_SYST; diagnosis group'),
 transformation_weight=c('A=0; B=1.263193; C=1.780242; D=1.514401; diagnosis-specific additions documented in methods','0.597904 if dependent; 0 if independent','1.576185 if ventilation; 0 if not','0.107441 * max(20 - BMI, 0)','A: 0.557670 * max(PASP - 40, 0)/10; B/C/D: 0.123048 * max(PASP, 20)/10'))
weights<-merge(weights,medians,by='component',sort=FALSE)
weights<-merge(weights,missing[era=='All',.(component,N,missing_n,missing_percent)],by='component',sort=FALSE)
weights[,missing_percent:=round(missing_percent,2)]
fwrite(weights,file.path(dest,'table1_component_definitions.csv'))
scores<-fread(file.path(dest,'table1_composite_summary.csv'))
scores[,`Median (IQR)`:=sprintf('%.2f (%.2f-%.2f)',median,q25,q75)]
updates<-fread(file.path(dest,'score_update_audit.csv'))
raw<-readRDS(file.path(revision_dir,'cache','lung_las_raw.rds'));b<-baseline_data()[WL_ORG=='LU']
raw<-merge(b[,.(PERS_ID,PX_ID)],raw,by=c('PERS_ID','PX_ID'),all.x=TRUE)
raw[,diagnosis_group:=las_diagnosis_group(CAN_DGN,CAN_PULM_ART_MEAN)]
raw[CAN_DGN %in% 1605 & is.na(CAN_PULM_ART_MEAN),diagnosis_group:=NA_character_]
raw[,functional_dependence:=fcase(CAN_FUNCTN_STAT %in% c(1,2070,2080,2090,2100,4070,4080,4090,4100),'Independent',CAN_FUNCTN_STAT %in% c(2,3,2010,2020,2030,2040,2050,2060,4010,4020,4030,4040,4050,4060),'Dependent',default='Missing')]
cats<-rbindlist(lapply(c('diagnosis_group','functional_dependence','CAN_VENTILATOR'),function(v){a<-raw[,.(n=.N),by=v];setnames(a,v,'category');a[,`:=`(variable=v,percent=round(100*n/nrow(raw),2))];a}))
fwrite(cats,file.path(dest,'table1_categorical_inputs.csv'))
continuous<-rbindlist(lapply(c('CAN_BMI','CAN_PULM_ART_SYST','CAN_PULM_ART_MEAN'),function(v){z<-as.numeric(raw[[v]]);z[!is.finite(z)|z<=0]<-NA_real_;data.table(variable=v,N=length(z),missing_n=sum(is.na(z)),missing_percent=round(100*mean(is.na(z)),2),median=median(z,na.rm=TRUE),q25=quantile(z,.25,na.rm=TRUE),q75=quantile(z,.75,na.rm=TRUE))}))
fwrite(continuous,file.path(dest,'table1_continuous_inputs.csv'))
esc<-function(x){x<-gsub('&','&amp;',as.character(x),fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE);gsub('>','&gt;',x,fixed=TRUE)}
htmltable<-function(x)paste0('<table><thead><tr>',paste0('<th>',esc(names(x)),'</th>',collapse=''),'</tr></thead><tbody>',paste(apply(as.data.frame(x),1,function(z)paste0('<tr>',paste0('<td>',esc(z),'</td>',collapse=''),'</tr>')),collapse=''),'</tbody></table>')
display<-r[,.(Model=model,Subgroup=subgroup,Pollutant,N,Events=events,`HR (95% CI)`=HR_CI,P,Status=status)]
tests[,P:=ifelse(p_value<.001,'<.001',sprintf('%.3f',p_value))]
intro<-paste0('<h1>Lung severity composite sensitivity</h1><p>The same study-defined clinical composite was used before and after March 9, 2023. It is <strong>neither LAS nor CAS</strong>; no organ-offer scores enter these models. Primary analyses and the original Table 1 are unchanged.</p>',
'<p>The composite sums five available domains using fixed weights from the 2021 LAS waiting-list mortality equation: diagnosis, functional dependence, ventilation, low BMI, and pulmonary artery systolic pressure. Age is adjusted separately. Cardiac index belongs to the final post-transplant component and was not added. Oxygen, PCO2, walk distance, bilirubin and creatinine were excluded for inadequate baseline availability or timing. Missing contributions use fixed observed pre-transition medians, with separate domain-missingness indicators in the models. Weights were not fitted to our outcomes.</p>',
'<p>Later dated registration forms update observed components; remaining values carry forward. Updates take effect the next day. Dense clinical histories are unavailable, so most candidates do not have repeated measurements. Higher composite values indicate higher modeled clinical risk under the borrowed weights, not a calibrated mortality probability or allocation score. Pediatric interpretation is exploratory.</p>',
'<p>Pooled models adjust for age, sex, race, ACS vulnerability, the composite and missingness indicators, stratified by first listing center. Listing year is not added. Exposure remains monthly time-varying: PM2.5 per 5 ug/m3, NO2 and O3 per 10 ppb. PM2.5/O3 end in 2024 and NO2 in 2025. Model-based standard errors and Efron ties match the original models.</p>')
parts<-paste0('<h2>Single- and multipollutant models</h2>',htmltable(display[Model %in% c('Single pollutant','PM2.5 + NO2','PM2.5 + NO2 + O3')]),
'<h2>Pre-/post-transition comparisons</h2><p>Listing-era comparisons group candidates by initial listing date; pre-CAS listings can have post-CAS follow-up. Calendar-era comparisons split follow-up on March 9, 2023, so candidates can contribute nonoverlapping time to both periods. The cutoff marks CAS implementation. The main era estimates below come from interaction models with center-by-era strata and common nuisance coefficients. They do not test a causal policy effect.</p>',htmltable(era_effects[,.(Pollutant=pollutant,Definition=era_definition,Era=era,`HR (95% CI)`,Interaction_P=interaction_P,Status=status)]),htmltable(tests),'<h3>Separate era fits (diagnostic)</h3><p>Separate fits allow nuisance effects to vary and can be unstable. Convergence warnings and failed fits below should not be interpreted as reliable estimates.</p>',htmltable(display[Model %in% c('listing_era','calendar_era')]),
'<h2>Table 1 addendum: score and components</h2>',htmltable(scores[,.(Era=era,N,`Median (IQR)`)]),htmltable(weights),htmltable(cats),htmltable(continuous),
'<h2>Availability and updates</h2>',htmltable(missing),htmltable(updates),htmltable(coverage),
'<h2>Other sensitivity and subgroup models</h2>',htmltable(display[!Model %in% c('Single pollutant','PM2.5 + NO2','PM2.5 + NO2 + O3','listing_era','calendar_era')]),
'<h2>Validation</h2><p>Original-score models fitted to the newly split intervals reproduced the original hazard ratios to within 0.000001 with identical candidate/event counts. Checks confirmed preserved person-time and events, contiguous intervals and no future score use.</p>',htmltable(rbindlist(validation)),
'<p>Full coefficients, formulas, warnings, imputation constants and component missingness are provided in companion CSV files. Exact weights and limitations are documented in docs/lung_severity_composite_sensitivity.md. LAS, lung allocation score; CAS, composite allocation score; BMI, body mass index; PASP, pulmonary artery systolic pressure; ACS, American Community Survey; HR, hazard ratio; CI, confidence interval; IQR, interquartile range.</p>')
out<-paste0('<!doctype html><html><head><meta charset="utf-8"><title>Lung severity composite sensitivity</title><style>body{font-family:Arial,sans-serif;color:#111;max-width:1400px;margin:35px auto;padding:0 24px;line-height:1.5}table{width:100%;border-collapse:collapse;margin:18px 0 32px;font-size:13px}th,td{text-align:left;padding:8px 10px;border-bottom:1px solid #ccc;vertical-align:top}th{border-top:2px solid black;border-bottom:2px solid black}h1{font-size:26px}h2{font-size:21px}tr{break-inside:avoid}thead{display:table-header-group}</style></head><body>',intro,parts,'</body></html>')
out<-gsub('PM2.5','PM<sub>2.5</sub>',out,fixed=TRUE);out<-gsub('NO2','NO<sub>2</sub>',out,fixed=TRUE);out<-gsub('O3','O<sub>3</sub>',out,fixed=TRUE)
writeLines(out,file.path(dest,'lung_composite_supplement.html'))
log_revision('Composite report and validation complete')
