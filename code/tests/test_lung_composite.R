source('code/revision_2026/lung_severity_composite.R')
d<-data.frame(CAN_DGN=rep(1607,4),CAN_PULM_ART_MEAN=25,CAN_FUNCTN_STAT=c(2070,2060,2070,NA),CAN_PULM_ART_SYST=c(30,30,50,NA),CAN_BMI=c(25,25,18,NA),ventilation=c(0,1,0,NA))
z<-lung_composite_components(d)
stopifnot(abs(z$functional_dependence[2]-.59790409246653)<1e-12,
 abs(z$ventilation[2]-1.57618530736936)<1e-12,
 abs(z$low_BMI[3]-.10744133677215*2)<1e-12,
 abs(z$pulmonary_pressure[3]-.55767046368853)<1e-12,
 z$diagnosis[1]==0,is.na(z$functional_dependence[4]))
reference<-vapply(z,median,na.rm=TRUE,numeric(1))
r<-impute_lung_components(z,reference)
stopifnot(all(is.finite(r$lung_composite)),r$missing_low_BMI[4]==1,r$missing_ventilation[4]==1,
 abs(r$lung_composite[2]-(.59790409246653+1.57618530736936))<1e-12)
d$CAN_DGN[1]<-999;d$CAN_DGN[2]<-1605;d$CAN_PULM_ART_MEAN[2]<-NA_real_
z<-lung_composite_components(d)
stopifnot(is.na(z$diagnosis[1]),is.na(z$diagnosis[2]),is.na(z$pulmonary_pressure[2]))
cat('Composite transformations, fixed weights, missingness and diagnosis tests passed.\n')
