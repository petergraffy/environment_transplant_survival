source("code/revision_2026/common.R")
suppressPackageStartupMessages(library(patchwork))
dat <- baseline_data()
dest <- file.path(revision_dir,"aj_with_risk_tables"); dir.create(dest,showWarnings=FALSE)
curves <- read_csv("output/prior_year_pollution_quartile_aalen_johansen_cif/prior_year_pollution_quartile_aalen_johansen_curve_data.csv",show_col_types=FALSE)
quartiles <- c("Q1 lowest","Q2","Q3","Q4 highest")
colors <- setNames(c("#2166AC","#67A9CF","#F4A582","#B2182B"),quartiles)
labels <- c("1st quartile (lowest)","2nd quartile","3rd quartile","4th quartile (highest)")
titles <- c("PM2.5"="PM[2.5]","NO2"="NO[2]","O3"="O[3]")
raws <- c("PM2.5"="pm25_prior_ug_m3","NO2"="no2_prior_ppb","O3"="o3_prior_ppb")
ticks <- function(h) switch(as.character(h),`1`=c(0,.5,1),`3`=0:3,`5`=0:5,`10`=seq(0,10,2))
counts <- list(); zeros <- list(); reconciliation <- list()
cohort_reference <- read_csv(file.path(revision_dir,"registration_spell_summary.csv"),show_col_types=FALSE)
original_aj <- read_csv("output/prior_year_pollution_quartile_aalen_johansen_cif/prior_year_pollution_quartile_event_summary.csv",show_col_types=FALSE)
for(p in names(raws)) {
  col <- raws[[p]]
  audit <- copy(dat[,c("WL_ORG","followup_days","adverse_event","transplant_or_improvement",col),with=FALSE])
  audit[,exclusion:=fcase(!is.finite(followup_days) | followup_days<=0,"Nonpositive or missing follow-up",
    is.na(adverse_event) | is.na(transplant_or_improvement),"Missing outcome status",
    !is.finite(get(col)),"Unavailable prelisting exposure",default="Included")]
  reconciliation[[p]] <- audit[,.(cohort_n=.N,
    excluded_followup=sum(exclusion=="Nonpositive or missing follow-up"),
    excluded_outcome=sum(exclusion=="Missing outcome status"),
    excluded_exposure=sum(exclusion=="Unavailable prelisting exposure"),
    aj_n=sum(exclusion=="Included")),by=WL_ORG][,pollutant:=p]
  d <- copy(dat[is.finite(get(col)) & followup_days>0 & !is.na(adverse_event) & !is.na(transplant_or_improvement)])
  br <- quantile(d[[col]],seq(0,1,.25),type=7); br[c(1,5)]<-c(-Inf,Inf)
  d[,quartile:=cut(get(col),br,labels=quartiles,include.lowest=TRUE)]
  for(org in names(organ_names)) for(q in quartiles) {
    x <- d[WL_ORG==org & quartile==q,followup_days/365.25]
    zeros[[length(zeros)+1L]] <- data.frame(organ=organ_names[[org]],pollutant=p,quartile=q,time=0,cif_adverse=0)
    for(h in c(1,3,5,10)) {
      tt<-ticks(h)
      counts[[length(counts)+1L]]<-data.frame(organ=organ_names[[org]],pollutant=p,quartile=q,horizon=h,time=tt,n_risk=risk_counts(x,tt))
    }
  }
}
risks <- bind_rows(counts)
stopifnot(all(risks$n_risk>=0),!anyDuplicated(risks[c("organ","pollutant","quartile","horizon","time")]))
stopifnot(all(vapply(split(risks,interaction(risks$organ,risks$pollutant,risks$quartile,risks$horizon)),function(x) all(diff(x$n_risk[order(x$time)])<=0),logical(1))))
write_csv(risks,file.path(dest,"numbers_at_risk.csv"))
zero <- risks %>% filter(horizon==1,time==0) %>% select(organ,pollutant,quartile,n_risk)
matched <- full_join(zero,original_aj,by=c("organ","pollutant","quartile"))
stopifnot(nrow(matched)==48L,all(matched$n_risk==matched$n))
audit <- bind_rows(reconciliation) %>% mutate(organ=recode(WL_ORG,!!!organ_names)) %>%
  left_join(cohort_reference %>% select(WL_ORG,reported_cohort_n=people),by="WL_ORG") %>%
  left_join(zero %>% group_by(organ,pollutant) %>% summarise(time_zero_total=sum(n_risk),.groups="drop"),by=c("organ","pollutant"))
stopifnot(all(audit$cohort_n==audit$reported_cohort_n),all(audit$aj_n==audit$time_zero_total),
  all(audit$cohort_n==audit$aj_n+audit$excluded_followup+audit$excluded_outcome+audit$excluded_exposure))
write_csv(audit,file.path(dest,"cohort_time_zero_reconciliation.csv"))
curves <- bind_rows(curves %>% select(organ,pollutant,quartile,time,cif_adverse),bind_rows(zeros)) %>% arrange(organ,pollutant,quartile,time)
common_theme <- theme_classic(base_size=24)+theme(panel.border=element_rect(fill=NA,color="black"),
  strip.background=element_rect(fill="white",color="black"),strip.text=element_text(face="bold",size=24),
  strip.text.y=element_text(angle=270),panel.spacing.x=unit(24,"pt"),
  axis.text=element_text(color="black"),legend.position="bottom",legend.title=element_blank())
# Recolor rendered axis labels by text so factor ordering cannot reverse the colors.
color_risk_labels <- function(plot) {
  palette <- setNames(unname(colors[quartiles]),c("1st","2nd","3rd","4th"))
  n_labels <- 0L
  recolor <- function(g) {
    if(inherits(g,"text") && is.character(g$label) && length(g$label)>0L &&
       all(g$label %in% names(palette))) {
      g$gp$col <- unname(palette[g$label])
      n_labels <<- n_labels + length(g$label)
    }
    if(length(g$grobs)) g$grobs <- lapply(g$grobs,recolor)
    if(length(g$children)) for(i in seq_along(g$children)) g$children[[i]] <- recolor(g$children[[i]])
    g
  }
  result <- recolor(patchwork::patchworkGrob(plot))
  stopifnot(n_labels==16L)
  result
}
for(h in c(1,3,5,10)) for(set in c("pm25_no2","o3")) {
  ps<-if(set=="pm25_no2") c("PM2.5","NO2") else "O3"
  rows<-list()
  for(org in c("Kidney","Liver","Heart","Lung")) {
    cd<-curves %>% filter(organ==.env$org,pollutant %in% ps,time<=h) %>% mutate(pollutant=factor(pollutant,levels=ps),quartile=factor(quartile,levels=quartiles))
    rd<-risks %>% filter(organ==.env$org,pollutant %in% ps,horizon==h) %>% mutate(pollutant=factor(pollutant,levels=ps),quartile=factor(quartile,levels=rev(quartiles)))
    a<-ggplot(cd,aes(time,cif_adverse,color=quartile))+geom_step(linewidth=1)+
      facet_grid(organ~pollutant,labeller=labeller(pollutant=as_labeller(titles,label_parsed)))+
      scale_color_manual(values=colors,breaks=quartiles,labels=labels)+
      scale_x_continuous(limits=c(0,h),breaks=ticks(h),expand=expansion(mult=c(.09,.09)))+
      scale_y_continuous(labels=scales::label_percent(accuracy=1),expand=expansion(mult=c(.02,.06)))+
      labs(x="Years after listing",y="Cumulative incidence")+common_theme
    b<-ggplot(rd,aes(time,quartile,label=scales::comma(n_risk)))+geom_text(size=5.4)+
      facet_grid(organ~pollutant,labeller=labeller(pollutant=as_labeller(titles,label_parsed)))+
      scale_x_continuous(limits=c(0,h),breaks=ticks(h),expand=expansion(mult=c(.09,.09)))+
      scale_y_discrete(labels=c("Q4 highest"="4th","Q3"="3rd","Q2"="2nd","Q1 lowest"="1st"))+
      labs(x=NULL,y="Number at risk")+common_theme+
      theme(strip.text.x=element_blank(),strip.text.y=element_text(color="transparent",angle=270,size=24,face="bold"),strip.background=element_blank(),panel.border=element_blank(),
            axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=18),
            plot.margin=margin(5,5,15,5))
    rows[[org]]<-(a/b)+plot_layout(heights=c(3.4,1.25))
  }
  plot<-wrap_plots(rows,ncol=1,guides="collect") & theme(legend.position="bottom")
  plot<-color_risk_labels(plot)
  width<-if(set=="pm25_no2") 19 else 11.5
  stem<-file.path(dest,paste0("aj_",set,"_",h,"yr_with_risk"))
  ggsave(paste0(stem,".png"),plot,width=width,height=33,dpi=250,limitsize=FALSE,bg="white")
  ggsave(paste0(stem,".pdf"),plot,width=width,height=33,limitsize=FALSE,bg="white")
}
writeLines("Numbers at risk count candidates still under observation and free of either adverse or competing events immediately before each time point. Quartile thresholds match the national, pollutant-specific prelisting exposure distribution. Counts beyond the final observed follow-up are zero.",file.path(dest,"risk_table_note.txt"))
