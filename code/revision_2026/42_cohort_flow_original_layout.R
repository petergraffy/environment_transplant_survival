source("code/revision_2026/common.R")
source("code/saf_paths.R")
library(haven)
setDTthreads(2)
out <- "output/figures/cohort_flow_diagram"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
cutoff <- as.Date("2023-03-08")
full <- baseline_data()
lung <- readRDS(file.path(revision_dir,"cache","las_extended_followup","prepared.rds"))$cohort
cohort <- rbindlist(list(full[WL_ORG != "LU"],lung),fill=TRUE)
components <- readRDS(file.path(revision_dir,"cache","primary_outcome_components.rds"))
cohort[, component := NULL]
cohort <- merge(cohort,components[,.(PERS_ID,WL_ORG,component)],by=c("PERS_ID","WL_ORG"),all.x=TRUE)
stopifnot(nrow(cohort)==1007801L,!anyNA(cohort$component))

# Define the eligible registration universe up front, not as a later lung exclusion.
# Later registrations belonging to a retained lung candidate remain part of deduplication.
paths <- get_saf_paths(release="q1_2026")
cols <- c("PERS_ID","PX_ID","WL_ORG","CAN_LISTING_DT","CAN_SOURCE","CAN_REM_DT","CAN_DEATH_DT","CAN_ENDWLFU")
raw <- rbindlist(lapply(c("cand_kipa.sas7bdat","cand_liin.sas7bdat","cand_thor.sas7bdat"),function(f)
  as.data.table(read_sas(file.path(paths$pubsaf_dir,f),col_select=all_of(cols)))),fill=TRUE)
zip <- as.data.table(read_sas(paths$canzip_file,col_select=all_of(c("PERS_ID","PX_ID","CAN_PERM_ZIP"))))
zip[, candidate_zip := stringr::str_extract(as.character(CAN_PERM_ZIP),"[0-9]+")]
zip <- unique(zip[,.(PERS_ID,PX_ID,candidate_zip)],by=c("PERS_ID","PX_ID"))
raw <- merge(raw,zip,by=c("PERS_ID","PX_ID"),all.x=TRUE)
raw[, index_date := CAN_LISTING_DT]
raw[, end_date := fcoalesce(as.Date(CAN_REM_DT),as.Date(CAN_DEATH_DT),as.Date(CAN_ENDWLFU),as.Date("2025-12-31"))]
raw[, end_date := pmin(end_date,as.Date("2025-12-31"))]
raw <- raw[WL_ORG %in% names(organ_names) & CAN_SOURCE %in% c("A","R") &
  !is.na(index_date) & index_date >= as.Date("2005-01-01") & index_date <= as.Date("2025-12-31") & end_date >= index_date]
missing_only_lung <- raw[WL_ORG=="LU",.(has_zip=any(!is.na(candidate_zip)),first=min(index_date)),by=PERS_ID][!has_zip & first<=cutoff,PERS_ID]
raw <- raw[WL_ORG!="LU" | PERS_ID %in% c(lung$PERS_ID,missing_only_lung)]
counts <- raw[!is.na(candidate_zip),.(registrations=.N),by=WL_ORG]
counts <- merge(counts,cohort[,.(candidates=.N),by=WL_ORG],by="WL_ORG")
counts[, duplicates := registrations-candidates]
counts[, organ := organ_names[WL_ORG]]
counts <- counts[match(c("HR","KI","LI","LU"),WL_ORG)]
n0 <- nrow(raw); nm <- sum(is.na(raw$candidate_zip)); nz <- n0-nm; nd <- sum(counts$duplicates)
stopifnot(nz-nd==nrow(cohort),all(counts$duplicates>=0),
  uniqueN(raw[!is.na(candidate_zip) & WL_ORG=="LU",PERS_ID])==nrow(lung))
summary <- cohort[,.(candidates=.N,adverse_events=sum(adverse_event),transplant=sum(component=="transplant"),
  median_waitlist_days=median(followup_days)),by=WL_ORG]
summary <- summary[match(c("HR","KI","LI","LU"),WL_ORG)]
summary[, organ_label := organ_names[WL_ORG]]
fmt <- function(x) format(x,big.mark=",",scientific=FALSE,trim=TRUE)
main <- data.table(id=c("eligible","zcta","candidates"),x=5.4,y=c(11,8.75,6.35),w=7.7,h=c(1.25,1.15,1.15),
  label=c(paste0("SRTR waitlist registrations meeting cohort eligibility\nHeart, kidney, liver, or lung\nN = ",fmt(n0)),
    paste0("Included registrations with residential ZCTA\nN = ",fmt(nz)),
    paste0("Deduplicated candidates\nN = ",fmt(nrow(cohort)))))
org <- data.table(id=paste0("organ_",summary$WL_ORG),x=c(1.25,4,6.8,9.55),y=4.35,w=2.35,h=1.1,
  label=paste0(summary$organ_label,"\nCandidates: ",fmt(summary$candidates)))
ends <- data.table(id=paste0("outcome_",summary$WL_ORG),x=org$x,y=1.45,w=2.35,h=2.2,
  label=paste0("Adverse events:\n",fmt(summary$adverse_events),"\nTransplants:\n",fmt(summary$transplant),
    "\nMedian waitlist:\n",fmt(summary$median_waitlist_days)," days"))
exc <- data.table(id=c("missing_zcta","duplicates"),x=10.95,y=c(9.9,7.55),w=3,h=c(1.15,2.25),
  label=c(paste0("Excluded\nMissing residential ZCTA\nn = ",fmt(nm)),
    paste0("Excluded\nDuplicate same-candidate,\nsame-organ registration rows\nn = ",fmt(nd),"\n",
      paste0(counts$organ,": ",fmt(counts$duplicates),collapse="\n"))))
boxes <- rbindlist(list(main,org,ends))
segments <- rbindlist(list(data.table(x=5.4,xend=5.4,y=c(10.375,8.175),yend=c(9.325,6.925)),
  data.table(x=5.4,xend=9.45,y=exc$y,yend=exc$y),
  data.table(x=5.4,xend=org$x,y=5.775,yend=4.9),
  data.table(x=org$x,xend=org$x,y=3.8,yend=2.55)))
p <- ggplot()+geom_segment(data=segments,aes(x=x,y=y,xend=xend,yend=yend),linewidth=.5,color="black")+
  geom_rect(data=boxes,aes(xmin=x-w/2,xmax=x+w/2,ymin=y-h/2,ymax=y+h/2),fill="white",color="black",linewidth=.65)+
  geom_rect(data=exc,aes(xmin=x-w/2,xmax=x+w/2,ymin=y-h/2,ymax=y+h/2),fill="#F2F2F2",color="black",linewidth=.65)+
  geom_text(data=boxes,aes(x,y,label=label),size=3.95,lineheight=.92)+
  geom_text(data=exc,aes(x,y,label=label),size=3.7,lineheight=.92)+
  coord_cartesian(xlim=c(0,12.8),ylim=c(.15,11.8),clip="off")+theme_void()+theme(plot.margin=margin(18,18,18,18))
ggsave(file.path(out,"cohort_flow_diagram.png"),p,width=13,height=11,dpi=600,bg="white")
ggsave(file.path(out,"cohort_flow_diagram.pdf"),p,width=13,height=11,device=cairo_pdf,bg="white")
fwrite(rbindlist(list(boxes,exc)),file.path(out,"cohort_flow_diagram_labels.csv"))
fwrite(summary,file.path(out,"cohort_flow_organ_outcomes.csv"))
fwrite(counts,file.path(out,"cohort_flow_registration_breakdown.csv"))
fwrite(data.table(step=c("Eligible registrations before ZCTA exclusion","Excluded: missing residential ZCTA",
  "Included: residential ZCTA available","Duplicate registrations removed","Deduplicated candidates"),n=c(n0,nm,nz,nd,nrow(cohort))),file.path(out,"cohort_flow_accounting.csv"))
writeLines(c("Eligibility is applied before the first box: heart/kidney/liver entry 2005-2025; lung entry through March 8, 2023.",
  "Later same-organ registrations for included lung candidates remain included in the duplicate registration accounting.",
  "Follow-up uses the extended-follow-up cohort, with outcomes through observed follow-up up to December 31, 2025.",
  "Counts are candidate-organ observations, not unique persons across organs. No separate lung-exclusion box is displayed."),file.path(out,"cohort_flow_note.txt"))
print(summary); print(counts)
