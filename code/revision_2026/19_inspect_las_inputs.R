source('code/revision_2026/common.R')
source('code/saf_paths.R')
p <- get_saf_paths(release='q1_2026')
f <- list.files(p$pubsaf_dir,pattern='^cand_thor[.]sas7bdat$',full.names=TRUE,ignore.case=TRUE)
h <- haven::read_sas(f,n_max=0)
fields <- unique(c('PERS_ID','PX_ID','WL_ORG','CAN_LISTING_DT',grep('CAN_(DGN|AGE|BMI|HGT|WGT|FUNCTN|ON_VENT|VENTILATOR|MED_COND|PCO2|AT_REST_O2|SIX_MIN|PULM_ART|CARDIAC_OUTPUT|MOST_RECENT_CREAT)',names(h),value=TRUE)))
d <- as.data.table(haven::read_sas(f,col_select=all_of(fields)))[WL_ORG=='LU']
dir.create(file.path(revision_dir,'lung_las_cas_sensitivity'),showWarnings=FALSE)
saveRDS(d,file.path(revision_dir,'cache','lung_las_raw.rds'))
print(names(d)); print(d[,.N,by=CAN_DGN][order(CAN_DGN)])
for(v in intersect(c('CAN_FUNCTN_STAT','CAN_MED_COND','CAN_ON_VENTILATOR','CAN_VENTILATOR'),names(d))) {cat('\n',v,'\n');print(d[,.N,by=v])}
b <- baseline_data(); cat('\nBASELINE FIELDS\n');print(names(b));cat('\nDATES\n');print(sapply(b[,intersect(names(b),c('index_date','observed_end_date')),with=FALSE],class))
