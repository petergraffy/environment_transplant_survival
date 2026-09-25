source("code/revision_2026/common.R")
base <- baseline_data()
cache <- file.path(revision_dir, "cache")
registration_file <- file.path(cache, "registrations.rds")
if (!file.exists(registration_file)) {
  log_revision("Reconstructing registration-level cohort using the existing SAF preparation")
  env <- load_script_prefix("code/50_primary_waitlist_period_pollution_cox.R", "registration_cohort")
  cols <- c("PERS_ID","PX_ID","WL_ORG","waitlist_row_id","index_date","observed_end_date",
            "candidate_zip","age","sex","race","listing_center","listing_year","adverse_event",
            "transplant_or_improvement","other_exit")
  reg <- as.data.table(env$cohort)[, ..cols]
  stopifnot(!anyDuplicated(reg[, .(PERS_ID, PX_ID, WL_ORG)]))
  setorder(reg, PERS_ID, WL_ORG, index_date, observed_end_date, waitlist_row_id)
  reg[, spell := spell_ids(index_date, observed_end_date), by = .(PERS_ID, WL_ORG)]
  saveRDS(reg, registration_file, compress = FALSE)
  rm(env); gc()
} else reg <- readRDS(registration_file)

counts <- reg[, .(registrations=.N, spells=max(spell), centers=uniqueN(listing_center)), by=.(PERS_ID, WL_ORG)]
write_csv(counts[, .(people=.N, registrations=sum(registrations), spells=sum(spells),
                     people_with_multiple_spells=sum(spells>1), extra_spells=sum(spells-1),
                     people_with_multiple_registrations=sum(registrations>1)), by=WL_ORG],
          file.path(revision_dir,"registration_spell_summary.csv"))
stopifnot(nrow(counts) == nrow(base))

spell_file <- file.path(cache, "spells_prepared.rds")
if (!file.exists(spell_file)) {
  log_revision("Collapsing overlapping registrations, preserving separate spells")
  spells <- reg[, .SD[1], by=.(PERS_ID,WL_ORG,spell)]
  # Resolve endpoints outside the one-row table to avoid data.table column scoping.
  ends <- reg[, {
    pr <- fifelse(adverse_event==1L,3L,fifelse(transplant_or_improvement==1L,2L,fifelse(other_exit==1L,1L,0L)))
    ix <- which(pr>0L)
    k <- if(length(ix)) tail(ix[order(observed_end_date[ix],pr[ix],waitlist_row_id[ix])],1L) else which.max(observed_end_date)
    .(end=observed_end_date[k],bad=adverse_event[k],good=transplant_or_improvement[k],other=other_exit[k])
  }, by=.(PERS_ID,WL_ORG,spell)]
  spells[ends,on=.(PERS_ID,WL_ORG,spell), `:=`(observed_end_date=i.end,adverse_event=i.bad,transplant_or_improvement=i.good,other_exit=i.other)]
  spells[, `:=`(followup_days=as.numeric(observed_end_date-index_date),
                 listing_year_int=as.integer(format(index_date,"%Y")))]
  prep <- load_script_prefix("code/81_prior_year_pollution_cox_svi.R", "analysis_dat")
  svi <- prep$make_complete_acs_svi_proxy(prep$community_path)
  spells <- as.data.table(left_join(as.data.frame(spells),svi,by=c("candidate_zip"="zip","listing_year_int"="analysis_year")))
  exposure <- read_rolling_prior_pollution(spells)
  for(p in names(pollutants)) {
    spells <- as.data.table(left_join(as.data.frame(spells),exposure[[p]],by=c("candidate_zip"="zip","index_date"="index_date")))
  }
  spells[, `:=`(pm25_prior_5ug=pm25_prior_ug_m3/5,no2_prior_10ppb=no2_prior_ppb/10,o3_prior_10ppb=o3_prior_ppb/10)]
  saveRDS(spells,spell_file,compress=FALSE)
}

concurrency_file <- file.path(cache,"concurrency.rds")
if(!file.exists(concurrency_file)) {
  log_revision("Building exact concurrent-listing change intervals")
  keys <- base[,.(PERS_ID,WL_ORG,primary_id=waitlist_row_id,start=index_date,end=observed_end_date)]
  joined <- merge(reg,keys,by=c("PERS_ID","WL_ORG"),all=FALSE,sort=FALSE)
  ints <- joined[, concurrent_intervals(index_date,observed_end_date,start[1],end[1]),by=.(primary_id)]
  saveRDS(ints,concurrency_file,compress=FALSE)
  desc <- ints[,.(maximum_concurrency=max(concurrent_listings),entry_concurrency=concurrent_listings[1],
                  days_without_open_listing=sum((tstop-tstart)[concurrent_listings==0L])),by=primary_id]
  desc <- merge(desc,keys[,.(primary_id,WL_ORG)],by="primary_id")
  write_csv(desc[,.(people=.N,multiple_at_entry=sum(entry_concurrency>1),
                    ever_concurrent=sum(maximum_concurrency>1),
                    any_gap=sum(days_without_open_listing>0),gap_person_days=sum(days_without_open_listing)),by=WL_ORG],
            file.path(revision_dir,"concurrent_listing_summary.csv"))
}
log_revision("Registration sensitivity preparation complete")
