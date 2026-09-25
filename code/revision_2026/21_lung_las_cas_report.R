source('code/revision_2026/common.R')
source('code/revision_2026/las_2021.R')
setDTthreads(2)
dest<-file.path(revision_dir,'lung_las_cas_sensitivity')
cache<-file.path(revision_dir,'cache','lung_las_cas')
d<-readRDS(file.path(cache,'scored_monthly_intervals.rds'))
terms<-c(pm25='pm25_interval_5ug',no2='no2_interval_10ppb',o3='o3_interval_10ppb')
benchfile<-file.path(dest,'validation_original_score_models.csv')
bench<-if(file.exists(benchfile))list(fread(benchfile)) else list()
for(p in if(file.exists(benchfile))character() else names(terms))for(era in c(FALSE,TRUE)) {
  log_revision('Validation fit: ',p,'; era adjustment=',era)
  end<-if(p=='no2')as.Date('2025-12-31') else as.Date('2024-12-31')
  z<-d[start_date<=end]
  rhs<-c(terms[[p]],'age','sex','race','zcta_svi_proxy','organ_score',if(era)'post_cas','strata(listing_center)')
  need<-c(terms[[p]],'age','sex','race','zcta_svi_proxy','organ_score','listing_center')
  z<-z[complete.cases(z[,..need])]
  fit<-coxph(as.formula(paste('Surv(tstart,tstop,tv_adverse_event)~',paste(rhs,collapse='+'))),data=as.data.frame(z),ties='efron',robust=FALSE)
  q<-as.data.table(broom::tidy(fit,exponentiate=TRUE,conf.int=TRUE))[term==terms[[p]]]
  q[,`:=`(pollutant=p,era_adjustment=era,N=uniqueN(z$PERS_ID),events=sum(z$tv_adverse_event))]
  bench[[length(bench)+1L]]<-q
}
bench<-rbindlist(bench)
old<-fread('output/timevarying_pollution_cox_svi/timevarying_pollution_cox_svi_results.csv')[organ=='LU']
chk<-merge(bench[era_adjustment==FALSE],old,by='pollutant')
stopifnot(all(abs(chk$estimate-chk$hazard_ratio)<1e-6),all(chk$N==chk$candidate_organ_episodes),all(chk$events==chk$adverse_events))
fwrite(bench,file.path(dest,'validation_original_score_models.csv'))

res<-fread(file.path(dest,'pollutant_results.csv'));logs<-fread(file.path(dest,'model_log.csv'))
res[,pollutant:=names(terms)[match(term,terms)]]
res<-merge(res,unique(logs[,.(model,subgroup,formula)]),by=c('model','subgroup'),all.x=TRUE,allow.cartesian=TRUE)
# Formula differs across pollutants: retain only the matching exposure formula.
res<-res[mapply(function(f,t)grepl(t,f,fixed=TRUE),formula,term)]
res<-unique(res,by=c('model','subgroup','term'))
res[,HR_CI:=sprintf('%.2f (%.2f-%.2f)',HR,HR_low,HR_high)]
res[,p:=ifelse(p.value<.001,'<.001',sprintf('%.3f',p.value))]
res[,status:=vapply(seq_len(.N),function(i){x<-logs[model==res$model[i]&subgroup==res$subgroup[i]];x<-x[vapply(formula,function(f)grepl(res$term[i],f,fixed=TRUE),logical(1))];paste(unique(x$status),collapse='; ')},character(1))]
fwrite(res,file.path(dest,'supplement_results.csv'))
fwrite(res,file.path(dest,'supplement_results.tsv'),sep='\t')
simple<-res[model %in% c('Single pollutant','PM2.5 + NO2','PM2.5 + NO2 + O3'),.(Model=model,Pollutant=pollutant,Candidates=candidates,Events=events,`HR (95% CI)`=HR_CI,P=p,Status=status)]
t1<-fread(file.path(dest,'table1_lung_las_cas_inputs.csv'))
t1[variable %in% c('diagnosis','functional_status','mechanical_ventilation'),median_IQR:='Categorical; see distribution below']
t1[,missing_percent:=round(missing_percent,1)]
raw<-readRDS(file.path(revision_dir,'cache','lung_las_raw.rds'))
b<-baseline_data()[WL_ORG=='LU']
raw<-merge(b[,.(PERS_ID,PX_ID)],raw,by=c('PERS_ID','PX_ID'),all.x=TRUE)
raw[,group:=las_diagnosis_group(CAN_DGN,CAN_PULM_ART_MEAN)]
cat1<-rbindlist(list(raw[,.(n=.N),by=.(variable=rep('Diagnosis group',nrow(raw)),level=ifelse(is.na(group),'Unmapped',group))],raw[,.(n=.N),by=.(variable=rep('Functional independence',nrow(raw)),level=fcase(is.na(CAN_FUNCTN_STAT)|CAN_FUNCTN_STAT %in% c(996,998),'Missing',CAN_FUNCTN_STAT %in% c(1,2070:2100,4070:4100),'No assistance',default='Some/total assistance'))],raw[,.(n=.N),by=.(variable=rep('Mechanical ventilation',nrow(raw)),level=ifelse(CAN_VENTILATOR==1,'Yes','No'))]))
cat1[,percent:=round(100*n/nrow(raw),1)]
fwrite(cat1,file.path(dest,'table1_categorical_inputs.csv'))
fwrite(t1,file.path(dest,'table1_lung_las_cas_inputs.csv'))
coverage<-fread(file.path(dest,'score_coverage.csv'));updates<-fread(file.path(dest,'score_updates.csv'))
esc<-function(x){x<-gsub('&','&amp;',as.character(x),fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE);gsub('>','&gt;',x,fixed=TRUE)}
table_html<-function(x)paste0('<table><thead><tr>',paste0('<th>',esc(names(x)),'</th>',collapse=''),'</tr></thead><tbody>',paste(apply(as.data.frame(x),1,function(r)paste0('<tr>',paste0('<td>',esc(r),'</td>',collapse=''),'</tr>')),collapse=''),'</tbody></table>')
intro<-paste0('<h1>Lung LAS/CAS sensitivity analysis</h1><p><strong>Exploratory reconstruction only. Primary models and the original Table 1 are unchanged.</strong></p>',
'<p>The 2021 LAS formula was reconstructed from available candidate fields before March 9, 2023. Inputs were carried forward until a subsequent dated listing supplied information. LAS and post-transition offer-level total CAS entered as separate continuous terms per 10 points, with a transition indicator and separate score-unavailable indicators. Score updates after listing took effect the following day. CAS was the median of distinct same-day observed scores, never filled backward from a later offer.</p>',
'<p><strong>Important limitation:</strong> dated creatinine and bilirubin were unavailable for lung candidates. Official component-specific missing-value defaults, including creatinine 0.1 mg/dL for waiting-list survival and 40 mg/dL for adult post-transplant survival, were used. Most oxygen, walk-distance and PCO2 fields were also unavailable. The resulting baseline LAS median was 0.39 (IQR, 0.24-1.95); this is a default-dominated research estimate, not a credible recovery of clinical allocation scores. CAS includes donor/offer-dependent allocation components and is not solely medical severity. These models should not be represented as adjustment for fully observed longitudinal LAS/CAS.</p>',
'<p>Exposure was unchanged: daily-derived monthly PM2.5 and O3 through 2024, and monthly NO2 where available with annual fallback through 2025. Covariates were age, sex, race, ACS-derived vulnerability proxy and listing-center strata; listing year was not added. Age subgroup models omitted age; sex/race subgroup models used attained age and omitted the stratifying covariate. Efron ties and model-based standard errors were retained. No original score was included alongside LAS/CAS. The same candidates, pollutant coverage, events and follow-up rules were retained. Models with convergence warnings are explicitly flagged below.</p>')
sections<-paste0('<h2>Single- and multipollutant results</h2>',table_html(simple),'<h2>Table 1 addendum: lung score inputs</h2>',table_html(t1),'<h3>Categorical inputs</h3>',table_html(cat1),'<h2>Score coverage</h2>',table_html(coverage),table_html(updates),'<h2>All sensitivity and subgroup results</h2>',table_html(res[,.(Model=model,Subgroup=subgroup,Pollutant=pollutant,N=candidates,Events=events,`HR (95% CI)`=HR_CI,P=p,Status=status)]),'<h2>Validation and era-adjustment benchmark</h2><p>Refitting the original score model on the newly split intervals reproduced the original estimates (absolute HR difference &lt;0.000001) and identical candidate/event counts. The additional transition-indicator benchmark distinguishes a policy-era adjustment from the score substitution.</p>',table_html(bench))
writeLines(paste0('<!doctype html><html><head><meta charset="utf-8"><title>Lung LAS/CAS sensitivity</title><style>body{font-family:Arial,sans-serif;color:#111;max-width:1400px;margin:40px auto;padding:0 25px;font-size:15px;line-height:1.5}table{border-collapse:collapse;width:100%;margin:20px 0 32px;font-size:13px}th,td{text-align:left;padding:9px 10px;border-bottom:1px solid #ccc;vertical-align:top}th{border-top:2px solid #111;border-bottom:2px solid #111}h1{font-size:26px}h2{font-size:21px;margin-top:38px}tr{break-inside:avoid}thead{display:table-header-group}@media print{body{margin:0;padding:0}table{font-size:10px}}</style></head><body>',intro,sections,'<p>Abbreviations: LAS, lung allocation score; CAS, composite allocation score; WL, waiting list; TX, post-transplant; CO, cardiac output; BMI, body mass index; PCO2, partial pressure of carbon dioxide; PA, pulmonary artery; HR, hazard ratio; CI, confidence interval; IQR, interquartile range; ACS, American Community Survey.</p></body></html>'),file.path(dest,'lung_las_cas_supplement.html'))
log_revision('Validation passed and supplement tables exported.')
