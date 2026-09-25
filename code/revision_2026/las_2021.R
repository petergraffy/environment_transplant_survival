# Final pre-CAS (September 2021) OPTN LAS equation, for sensitivity use only.
# See data/reference/las_2021/README.md for provenance and limitations.
las_diagnosis_group <- function(code, pa_mean) {
  g <- rep(NA_character_, length(code))
  g[code %in% c(100,103,105,107,109,111,113,114,116,214,412,1603,1606,1607,1608,1611)] <- 'A'
  g[code %in% c(200,202,203,205,206,208,209,210,212,216,218,220,1500:1502,1517,1548,1549,1600,1601,1614,1615)] <- 'B'
  g[code %in% c(300,302,303,305,1602)] <- 'C'
  g[code %in% c(106,213,215,217,219,400:409,411,413:424,432,434,437,438,440,441,444,446:449,451,453,1518:1525,1550:1557,1599,1604,1609,1610,1612,1613)] <- 'D'
  # Sarcoidosis with unavailable mean pressure belongs to A under 2021 policy.
  g[code %in% 1605] <- ifelse(!is.na(pa_mean[code %in% 1605]) & pa_mean[code %in% 1605]>30,'D','A')
  g
}

las_auc <- function(lp, baseline) {
  vapply(lp, function(x) if(is.finite(x)) sum(baseline ^ exp(x)) else NA_real_, numeric(1))
}

calculate_las_2021 <- function(d) {
  n <- nrow(d)
  value <- function(v, default) {x<-as.numeric(d[[v]]); x[!is.finite(x)]<-default; x}
  age <- d$score_age
  g <- las_diagnosis_group(d$CAN_DGN,d$CAN_PULM_ART_MEAN)
  wl_group <- c(A=0,B=1.26319338239175,C=1.78024171092307,D=1.51440083414275)[g]
  tx_group <- c(A=0,B=.51341349576197,C=.23187885123342,D=.12527366545917)[g]
  bmi <- value('CAN_BMI',100)
  o2_wl <- value('CAN_AT_REST_O2',0); o2_tx<-value('CAN_AT_REST_O2',26.33)
  walk_wl<-value('CAN_SIX_MIN_WALK',4000);walk_tx<-value('CAN_SIX_MIN_WALK',0)
  pco2<-pmax(value('CAN_PCO2',40),40)
  pas<-pmax(value('CAN_PULM_ART_SYST',20),20)
  ci<-value('cardiac_index',3)
  vent_wl<-value('ventilation',0);vent_tx<-value('ventilation',1)
  # No temporally verified creatinine/bilirubin or serial PCO2 rise in these LU files.
  creat_wl<-value('dated_creatinine',.1);creat_tx<-value('dated_creatinine',40)
  bili<-pmax(value('dated_bilirubin',.7),.7)
  fs<-d$CAN_FUNCTN_STAT
  independent<-is.na(fs) | fs %in% c(996,998,1,2070:2100,4070:4100)
  bronch<-as.integer(d$CAN_DGN %in% 1608 & g=='A')
  fib<-as.integer(d$CAN_DGN %in% 1613 & g=='D')
  oblit<-as.integer(d$CAN_DGN %in% c(106,1612) & g=='D')
  sar_a<-as.integer(d$CAN_DGN %in% 1605 & g=='A')
  sar_d<-as.integer(d$CAN_DGN %in% 1605 & g=='D')
  wl <- .0281444188123287*age + .15572123729572*pmax(bili-1,0) +
    .10744133677215*pmax(20-bmi,0) + 1.57618530736936*vent_wl +
    .0996197163645*creat_wl*(age>=18) + wl_group + .40107198445555*bronch +
    .2088684500011*fib - .64590852776042*sar_d + 1.39885489102977*sar_a -
    .59790409246653*independent + ifelse(g=='B',.0340531822566417,.08232292818591)*o2_wl +
    .12639905519026*pco2/10 + ifelse(g=='A',.55767046368853*pmax(pas-40,0)/10,.1230478043299*pas/10) -
    .09937981549564*walk_wl/100
  tx <- .0208895939056676*pmax(age-45,0) + .25451764981323*creat_tx*(age>=18) +
    .1448727551614*(ci<2) + .33161555489537*vent_tx + tx_group + .12048575705296*bronch -
    .33402539276216*oblit + .43537371336129*sar_d + .98051166673574*sar_a +
    ifelse(g=='A',.0100383613234584,.0093694370076423)*o2_tx +
    .0001943695814883*pmax(1200-walk_tx,0)
  base<-read.csv('data/reference/las_2021/baseline_survival.csv')
  stopifnot(identical(base$day,0:364), all(diff(base[[2]])<=0), all(diff(base[[3]])<=0))
  wl_auc<-las_auc(wl,base[[2]]);tx_auc<-las_auc(tx,base[[3]])
  score<-100*(tx_auc-2*wl_auc+730)/1095
  score[is.na(age)|age<12]<-NA_real_
  data.frame(las=score,diagnosis_group=g,wl_lp=wl,tx_lp=tx,wl_auc,tx_auc)
}
