source("code/revision_2026/common.R")
suppressPackageStartupMessages(library(patchwork))
dest <- file.path(revision_dir,"svi_validation")
d <- read_csv(file.path(dest,"zcta_validation_values.csv"),show_col_types=FALSE)
stats <- read_csv(file.path(dest,"validation_results.csv"),show_col_types=FALSE)
long <- bind_rows(lapply(c("overall","socioeconomic"),function(metric) {
  d %>% filter(is.finite(zcta_svi_proxy),is.finite(.data[[metric]])) %>%
    transmute(metric,x=zcta_svi_proxy,y=.data[[metric]])
}))
stopifnot(nrow(long)==2L*32428L,all(long$x>=0 & long$x<=1),all(long$y>=0 & long$y<=1))
# Fixed 0.02-wide bins give both panels exactly the same spatial and color scale.
bins <- long %>% mutate(ix=pmin(49L,floor(x/.02)),iy=pmin(49L,floor(y/.02))) %>%
  count(metric,ix,iy,name="n") %>% mutate(x=(ix+.5)*.02,y=(iy+.5)*.02)
stopifnot(all((bins %>% group_by(metric) %>% summarise(n=sum(n)))$n==32428L))
upper <- ceiling(max(bins$n)/100)*100
breaks <- sort(unique(c(0,pretty(c(0,upper),n=5),upper)))
breaks <- breaks[breaks>=0 & breaks<=upper]
panels <- lapply(c("overall","socioeconomic"),function(metric) {
  s <- stats %>% filter(cdc_measure==metric)
  rho <- sprintf("Spearman rho = %.3f",s$spearman)
  ggplot(bins %>% filter(.data$metric==.env$metric),aes(x,y,fill=n))+
    geom_tile(width=.02,height=.02)+
    annotate("text",x=.035,y=.965,label=rho,hjust=0,vjust=1,size=5.2,color="black")+
    scale_fill_gradientn(colors=c("#eff5f8","#a6cddd","#4f96b8","#19577d","#082b49"),
      limits=c(0,upper),breaks=breaks,trans="sqrt",labels=scales::label_comma(),
      name="Number of ZCTAs per bin",
      guide=guide_colorbar(title.position="top",title.hjust=.5,
        barwidth=unit(16,"cm"),barheight=unit(.5,"cm"),
        ticks=TRUE,frame.colour="black",ticks.colour="black"))+
    scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.2),expand=expansion(mult=0))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2),expand=expansion(mult=0))+
    coord_fixed(clip="off")+
    labs(title=if(metric=="overall") "Overall SVI" else "Socioeconomic theme",
      x="ACS-derived vulnerability index",y="CDC SVI percentile rank")+
    theme_classic(base_size=17)+
    theme(panel.border=element_rect(color="black",fill=NA,linewidth=.6),
      axis.line=element_blank(),axis.text=element_text(color="black"),
      axis.title=element_text(size=17),axis.title.x=element_text(margin=margin(t=10)),
      axis.title.y=element_text(margin=margin(r=10)),
      plot.title=element_text(size=20,face="bold",hjust=0),
      legend.title=element_text(size=16,face="bold"),legend.text=element_text(size=14),
      plot.tag=element_text(size=20,face="bold"),plot.margin=margin(12,18,8,12))
})
fig <- wrap_plots(panels,nrow=1,guides="collect")+
  plot_annotation(tag_levels="A") & theme(legend.position="bottom")
stem <- file.path(dest,"svi_validation_AB")
ggsave(paste0(stem,".png"),fig,width=13.5,height=7.2,dpi=400,bg="white")
ggsave(paste0(stem,".pdf"),fig,width=13.5,height=7.2,bg="white")
writeLines(c("Supplemental Figure. Agreement of the ACS-derived neighborhood vulnerability index with CDC SVI in 2022.",
  "Panel A compares the index with overall CDC SVI; panel B compares it with the socioeconomic-status theme. Each panel includes 32,428 ZCTAs with valid values for both measures. Higher values indicate greater vulnerability. Shading indicates ZCTA counts within fixed 0.02-by-0.02 bins; both panels share a square-root-transformed color scale labeled in raw counts. Empty bins are white. Spearman rank correlations are displayed. Comparisons are unweighted by population or candidate count and use observed 2022 ACS data, not carried-forward or carried-backward vintages.",
  "ACS, American Community Survey; CDC, Centers for Disease Control and Prevention; SVI, Social Vulnerability Index; ZCTA, ZIP Code Tabulation Area."),
  paste0(stem,"_caption.txt"))
