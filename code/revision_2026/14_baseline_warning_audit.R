source("code/revision_2026/common.R")
audit_rows <- list()
audit_coxph <- function(...) {
  cl <- match.call(expand.dots=TRUE);cl[[1]] <- quote(survival::coxph)
  messages <- character()
  fit <- withCallingHandlers(eval(cl,envir=parent.frame()),warning=function(w) {
    messages <<- c(messages,conditionMessage(w));invokeRestart("muffleWarning")
  })
  for(msg in messages) {
    indices <- sub(".*variable +([0-9, ]+) *;.*","\\1",msg)
    idx <- suppressWarnings(as.integer(trimws(strsplit(indices,",",fixed=TRUE)[[1]])))
    terms <- if(all(!is.na(idx))) paste(names(coef(fit))[idx],collapse="; ") else NA_character_
    audit_rows[[length(audit_rows)+1L]] <<- data.frame(analysis=fit_context,organ=org,pollutant=p,message=msg,affected_terms=terms)
  }
  fit
}
coxph <- audit_coxph
for(script in c("03_baseline_diagnostics_eras.R","06_recent_finegray.R","09_recent_cause_specific.R")) {
  if(script=="06_recent_finegray.R") fit_context <- "recent_FineGray"
  if(script=="09_recent_cause_specific.R") fit_context <- "recent_cause_specific"
  statements <- as.list(parse(file.path("code/revision_2026",script)))
  for(statement in statements) {
    if(is.call(statement) && identical(statement[[1]],as.name("for"))) {
      dest <- file.path(revision_dir,"baseline_warning_audit",sub("[.]R$","",script))
      dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    }
    # Keep the comparator's final output in the audit directory as well.
    if(script=="09_recent_cause_specific.R" && is.call(statement) && identical(statement[[1]],as.name("write_csv"))) next
    eval(statement,envir=.GlobalEnv)
  }
}
write_csv(bind_rows(audit_rows),file.path(revision_dir,"baseline_warning_audit","captured_warnings.csv"))
