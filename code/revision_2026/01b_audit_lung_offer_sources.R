source("code/revision_2026/common.R")
source("code/saf_paths.R")
paths <- get_saf_paths(release="q1_2026")
files <- list.files(file.path(paths$saf_dir,"ptr_lu_2018_2025"),pattern="sas7bdat$",full.names=TRUE)
schema <- bind_rows(lapply(files,function(path) {
  d <- haven::read_sas(path,n_max=0)
  data.frame(file=basename(path),variable=names(d),label=vapply(d,function(x){z<-attr(x,"label");if(is.null(z)) "" else z},character(1)))
}))
write_csv(schema,file.path(revision_dir,"lung_offer_variable_schema.csv"))
print(as_tibble(schema %>% filter(grepl("score|LAS|CAS|date|match",paste(variable,label),ignore.case=TRUE))),n=120)
history <- haven::read_sas(file.path(paths$pubsaf_dir,"stathist_thor.sas7bdat"),
  col_select=any_of(c("WL_ORG","CANHX_STAT_CD","CANHX_BEGIN_DT","CANHX_END_DT")))
if("WL_ORG" %in% names(history)) history <- filter(history,WL_ORG=="LU")
write_csv(history %>% count(CANHX_STAT_CD,sort=TRUE),file.path(revision_dir,"lung_history_status_codes.csv"))
