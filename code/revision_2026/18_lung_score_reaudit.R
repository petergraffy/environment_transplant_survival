source("code/revision_2026/common.R")
source("code/saf_paths.R")
paths <- get_saf_paths(release="q1_2026")
dest <- file.path(revision_dir,"lung_score_reaudit")
dir.create(dest,showWarnings=FALSE)
files <- list.files(paths$saf_dir,pattern="[.]sas7bdat$",recursive=TRUE,full.names=TRUE)
files <- files[!grepl("/ptr_",gsub("\\\\","/",files)) | grepl("/ptr_lu_2018_2025/",gsub("\\\\","/",files))]
schema <- bind_rows(lapply(files,function(f) {
  log_revision("Header: ",basename(f))
  x <- haven::read_sas(f,n_max=0)
  tibble(file=substring(f,nchar(paths$saf_dir)+2L),variable=names(x),
    label=vapply(x,function(v) {s<-attr(v,"label");if(is.null(s)) "" else s},character(1)))
}))
write_csv(schema,file.path(dest,"fresh_schema.csv"))
hits <- schema %>% filter(grepl("score|wlauc|survival.*measure|urgency|(^|_)las($|_)|(^|_)cas($|_)",paste(variable,label),ignore.case=TRUE))
write_csv(hits,file.path(dest,"score_fields.csv"))
ptr <- files[grepl("ptr_lu_2018_2025/",gsub("\\\\","/",files))]
out <- list(); statuses <- list()
for(f in ptr) {
  log_revision("Score profile: ",basename(f))
  d <- as.data.table(haven::read_sas(f,col_select=c(PX_ID,MATCH_SUBMIT_DT,PTR_TOT_SCORE,PTR_STAT_CD)))
  d[,day:=as.Date(MATCH_SUBMIT_DT)]
  d[,era:=ifelse(day<as.Date("2023-03-09"),"Before CAS","CAS period")]
  for(e in unique(d$era)) {
    a <- d[era==e]; v <- a$PTR_TOT_SCORE[is.finite(a$PTR_TOT_SCORE)]
    statuses[[length(statuses)+1L]] <- a[, .N,by=PTR_STAT_CD][,`:=`(file=basename(f),era=e)]
    daily <- if(length(v)) a[is.finite(PTR_TOT_SCORE),.(distinct_scores=uniqueN(PTR_TOT_SCORE),range=max(PTR_TOT_SCORE)-min(PTR_TOT_SCORE)),by=.(PX_ID,day)] else data.table(distinct_scores=integer(),range=numeric())
    out[[length(out)+1L]] <- tibble(file=basename(f),era=e,rows=nrow(a),registrations=uniqueN(a$PX_ID),
      missing_score=sum(!is.finite(a$PTR_TOT_SCORE)),min=if(length(v)) min(v) else NA_real_,p25=if(length(v)) quantile(v,.25) else NA_real_,median=if(length(v)) median(v) else NA_real_,p75=if(length(v)) quantile(v,.75) else NA_real_,max=if(length(v)) max(v) else NA_real_,
      registration_days=nrow(daily),days_with_multiple_scores=sum(daily$distinct_scores>1),
      days_with_range_over_one=sum(daily$range>1))
  }
  rm(d);gc()
}
write_csv(bind_rows(out),file.path(dest,"match_score_profile.csv"))
write_csv(bind_rows(statuses),file.path(dest,"match_status_values.csv"))
histfile <- files[grepl("/stathist_thor.sas7bdat$",gsub("\\\\","/",files))]
h <- as.data.table(haven::read_sas(histfile))
write_csv(h[WL_ORG=="LU", .N,by=CANHX_STAT_CD][order(CANHX_STAT_CD)],file.path(dest,"lung_status_values.csv"))
