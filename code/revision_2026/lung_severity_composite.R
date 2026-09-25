source('code/revision_2026/las_2021.R')

# Reduced 2021 LAS waiting-list linear predictor, not LAS or CAS.
# Age remains a separate regression covariate; no post-transplant terms enter.
lung_composite_components <- function(d) {
  g<-las_diagnosis_group(d$CAN_DGN,d$CAN_PULM_ART_MEAN)
  g[d$CAN_DGN %in% 1605 & is.na(d$CAN_PULM_ART_MEAN)]<-NA_character_
  diagnosis<-unname(c(A=0,B=1.26319338239175,C=1.78024171092307,D=1.51440083414275)[g])
  diagnosis<-diagnosis + .40107198445555*(d$CAN_DGN %in% 1608 & g=='A') +
    .2088684500011*(d$CAN_DGN %in% 1613 & g=='D') -
    .64590852776042*(d$CAN_DGN %in% 1605 & g=='D') +
    1.39885489102977*(d$CAN_DGN %in% 1605 & g=='A')
  fs<-as.numeric(d$CAN_FUNCTN_STAT)
  independent<-fs %in% c(1,2070,2080,2090,2100,4070,4080,4090,4100)
  dependent<-fs %in% c(2,3,2010,2020,2030,2040,2050,2060,4010,4020,4030,4040,4050,4060)
  function_term<-ifelse(independent,0,ifelse(dependent,.59790409246653,NA_real_))
  pas<-as.numeric(d$CAN_PULM_ART_SYST)
  pas[!is.finite(pas)|pas<=0]<-NA_real_
  pressure<-ifelse(g=='A',.55767046368853*pmax(pas-40,0)/10,.1230478043299*pmax(pas,20)/10)
  bmi<-as.numeric(d$CAN_BMI);bmi[!is.finite(bmi)|bmi<=0]<-NA_real_
  vent<-as.numeric(d$ventilation);vent[!vent %in% c(0,1)]<-NA_real_
  data.frame(diagnosis=diagnosis,functional_dependence=function_term,
    ventilation=1.57618530736936*vent,low_BMI=.10744133677215*pmax(20-bmi,0),pulmonary_pressure=pressure)
}

impute_lung_components <- function(components, reference_medians) {
  stopifnot(identical(names(components),names(reference_medians)),all(is.finite(reference_medians)))
  missing<-as.data.frame(lapply(components,function(x)as.integer(!is.finite(x))))
  names(missing)<-paste0('missing_',names(components))
  for(v in names(components))components[[v]][!is.finite(components[[v]])]<-reference_medians[[v]]
  data.frame(lung_composite=rowSums(components),missing)
}
