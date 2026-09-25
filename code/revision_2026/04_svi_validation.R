source("code/revision_2026/common.R")
env <- load_script_prefix("code/81_prior_year_pollution_cox_svi.R","analysis_dat")
dest <- file.path(revision_dir,"svi_validation"); dir.create(dest,showWarnings=FALSE)
proxy <- env$make_complete_acs_svi_proxy(env$community_path) %>% filter(analysis_year==2022)
original <- read_csv(env$community_path,show_col_types=FALSE) %>% filter(analysis_year==2022)
stopifnot(all(original$community_year==2022),!anyDuplicated(original$zip))
components <- c("pct_poverty","pct_unemployed","pct_no_vehicle","pct_nonwhite","median_household_income","pct_bachelor_plus")
write_csv(tibble(component=components,n_zctas=nrow(original),missing=vapply(original[components],function(x) sum(is.na(x)),integer(1))),file.path(dest,"proxy_component_missingness.csv"))
write_csv(original %>% count(analysis_year,community_year),file.path(dest,"observed_acs_vintage.csv"))
write_csv(tibble(variable=names(original)),file.path(dest,"acs_available_fields.csv"))
svi <- read_csv("data/reference/cdc_svi/SVI_2022_US_ZCTA.csv",col_types=cols(FIPS=col_character()),show_col_types=FALSE) %>%
  transmute(zip=stringr::str_pad(FIPS,5,pad="0"),overall=if_else(RPL_THEMES<0,NA_real_,RPL_THEMES),
            socioeconomic=if_else(RPL_THEME1<0,NA_real_,RPL_THEME1),population=E_TOTPOP)
stopifnot(!anyDuplicated(proxy$zip),!anyDuplicated(svi$zip))
joined <- full_join(proxy,svi,by="zip")
write_csv(tibble(proxy_zctas=nrow(proxy),cdc_zctas=nrow(svi),matched=sum(!is.na(joined$analysis_year)&!is.na(joined$population)),
                 proxy_without_cdc=sum(!is.na(joined$analysis_year)&is.na(joined$population)),
                 cdc_without_proxy=sum(is.na(joined$analysis_year)&!is.na(joined$population))),file.path(dest,"coverage.csv"))
result <- list(); plots <- list()
for(metric in c("overall","socioeconomic")) {
  d <- joined %>% filter(is.finite(zcta_svi_proxy),is.finite(.data[[metric]]))
  # Rank both measures among the same matched ZCTAs; preserve ties.
  q <- function(x) pmin(4L,as.integer(ceiling(4*rank(x,ties.method="average")/length(x))))
  a <- q(d$zcta_svi_proxy); b <- q(d[[metric]])
  tab <- table(factor(a,levels=1:4),factor(b,levels=1:4)); pr <- tab/sum(tab)
  weights <- abs(outer(1:4,1:4,"-"))/3
  kappa <- 1-sum(weights*pr)/sum(weights*outer(rowSums(pr),colSums(pr)))
  rho <- cor(d$zcta_svi_proxy,d[[metric]],method="spearman")
  set.seed(20260923)
  boot <- replicate(500,{ix<-sample.int(nrow(d),replace=TRUE);cor(d$zcta_svi_proxy[ix],d[[metric]][ix],method="spearman")})
  ci <- quantile(boot,c(.025,.975))
  result[[metric]] <- tibble(cdc_measure=metric,n_zctas=nrow(d),spearman=rho,rho_low=ci[1],rho_high=ci[2],
    pearson=cor(d$zcta_svi_proxy,d[[metric]]),quartile_agreement=mean(a==b),linear_weighted_kappa=kappa,
    top_quartile_jaccard=sum(a==4&b==4)/sum(a==4|b==4),top_quartile_sensitivity=sum(a==4&b==4)/sum(b==4))
  write_csv(as.data.frame(tab),file.path(dest,paste0(metric,"_quartile_agreement.csv")))
  p <- ggplot(d,aes(x=zcta_svi_proxy,y=.data[[metric]]))+geom_bin_2d(bins=50)+
    scale_fill_gradient(low="#e2eef2",high="#17436b",trans="sqrt",name="ZCTAs")+
    labs(x="ACS vulnerability proxy",y=paste("CDC SVI",metric),title=paste0("2022: Spearman rho = ",sprintf("%.2f",rho)))+
    theme_classic(base_size=17)
  ggsave(file.path(dest,paste0("svi_",metric,".png")),p,width=8,height=6,dpi=300)
  ggsave(file.path(dest,paste0("svi_",metric,".pdf")),p,width=8,height=6)
}
write_csv(bind_rows(result),file.path(dest,"validation_results.csv"))
write_csv(joined,file.path(dest,"zcta_validation_values.csv"))
writeLines(c("CDC source: https://svi.cdc.gov/Documents/Data/2022/csv/zcta/SVI_2022_US_ZCTA.csv",
 "SHA256: e2007404c6e63aca40fcfab64893599de8e3cd79887d61463b36b7ca60f2791a",
 "One observation per ZCTA; 2022 ACS proxy and CDC national 2022 SVI. No patient-level data.",
 "Bootstrap confidence intervals resample ZCTAs (500 replicates); spatial dependence is not modeled.",
 "Quartiles are ranked within the common matched geography; agreement is convergent validity, not interchangeability."),file.path(dest,"methods.txt"))
