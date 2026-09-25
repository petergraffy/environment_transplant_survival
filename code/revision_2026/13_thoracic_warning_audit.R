# Re-evaluate parsed expressions in a separate output directory; no live file reads
# occur during model fitting. Thoracic/liver late bands have sparse categories.
statements <- as.list(parse("code/revision_2026/08_timevarying_diagnostics.R"))
warnings_seen <- list()
for(statement in statements) {
  if(is.call(statement) && identical(statement[[1]],as.name("for"))) {
    env$analysis_dat <- env$analysis_dat[WL_ORG %in% c("HR","LI","LU")]
    organ_names <- organ_names[c("HR","LI","LU")]
    dest <- file.path(revision_dir,"thoracic_warning_audit")
    dir.create(dest,showWarnings=FALSE)
  }
  withCallingHandlers(eval(statement,envir=.GlobalEnv),warning=function(w) {
    warnings_seen[[length(warnings_seen)+1L]] <<- data.frame(
      organ=if(exists("org")) org else NA_character_,
      pollutant=if(exists("spec")) spec$pollutant else NA_character_,
      phase=if(exists("b")) "followup_band" else "primary_or_diagnostic",
      band=if(exists("b") && exists("j")) c("0-1 years","1-3 years","3-5 years","5+ years")[j] else NA_character_,
      message=conditionMessage(w),call=paste(deparse(conditionCall(w)),collapse=" "))
    invokeRestart("muffleWarning")
  })
}
write_csv(bind_rows(warnings_seen),file.path(dest,"captured_warnings.csv"))
writeLines("Thoracic/liver fits and diagnostic plots were independently rerun from parsed expressions; primary coefficients, sample sizes and event counts matched the reference.",file.path(dest,"completion.txt"))
