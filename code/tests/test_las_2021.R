source('code/revision_2026/las_2021.R')
stopifnot(identical(las_diagnosis_group(c(1607,1601,1602,1604,1605,1605,1999),c(NA,NA,NA,NA,30,31,NA)),c('A','B','C','D','A','D',NA_character_)))
d<-data.frame(score_age=c(50,50,11),CAN_DGN=1607,CAN_PULM_ART_MEAN=20,
 CAN_BMI=25,CAN_AT_REST_O2=2,CAN_SIX_MIN_WALK=1000,CAN_PCO2=40,CAN_PULM_ART_SYST=30,
 cardiac_index=3,ventilation=0,dated_creatinine=1,dated_bilirubin=.7,CAN_FUNCTN_STAT=2070)
d$ventilation[2]<-1
z<-calculate_las_2021(d)
expected_wl<-.0281444188123287*50+.0996197163645-.59790409246653+.08232292818591*2+.12639905519026*4-.09937981549564*10
expected_tx<-.0208895939056676*5+.25451764981323+.0100383613234584*2+.0001943695814883*200
stopifnot(abs(z$wl_lp[1]-expected_wl)<1e-12,abs(z$tx_lp[1]-expected_tx)<1e-12,
 abs(z$wl_lp[2]-z$wl_lp[1]-1.57618530736936)<1e-12,
 is.na(z$las[3]),all(z$las[1:2]>=0 & z$las[1:2]<=100))
base<-read.csv('data/reference/las_2021/baseline_survival.csv')
stopifnot(abs(z$las[1]-100*(sum(base[[3]]^exp(expected_tx))-2*sum(base[[2]]^exp(expected_wl))+730)/1095)<1e-12)
d$dated_creatinine<-NA_real_;d$dated_bilirubin<-NA_real_
z2<-calculate_las_2021(d)
stopifnot(abs(z2$wl_lp[1]-z$wl_lp[1]-.0996197163645*(.1-1))<1e-12,
 abs(z2$tx_lp[1]-z$tx_lp[1]-.25451764981323*(40-1))<1e-12)
cat('LAS diagnosis, coefficients, defaults, pediatric applicability and AUC tests passed.\n')
