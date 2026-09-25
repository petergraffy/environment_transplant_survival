# Preserve the original 2020-2022 sensitivity; use a distinct output directory.
options(revision.recent_start="2015-01-01",revision.recent_end="2019-12-31",
        revision.recent_dir="finegray_2015_2019")
source("code/revision_2026/06_recent_finegray.R")
rm(dat);gc()
source("code/revision_2026/09_recent_cause_specific.R")
fg <- read_csv(file.path(dest,"finegray_results.csv"),show_col_types=FALSE)
cs <- read_csv(file.path(dest,"matched_recent_cause_specific.csv"),show_col_types=FALSE)
check <- inner_join(fg,cs,by=c("organ","pollutant"),suffix=c("_fg","_cs"))
stopifnot(nrow(check)==12L,all(check$n_fg==check$n_cs),all(check$events_fg==check$events_cs),
          all(is.finite(fg$hazard_ratio)),all(is.finite(cs$hazard_ratio)))
writeLines("PASS: all 12 Fine-Gray and matched cause-specific models have identical eligible candidate and adverse-event counts and finite pollutant estimates.",file.path(dest,"validation.txt"))
