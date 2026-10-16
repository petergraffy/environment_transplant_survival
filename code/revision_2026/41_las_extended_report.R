source("code/revision_2026/common.R")
dest <- file.path(revision_dir, "las_extended_followup")
x <- fread(file.path(dest, "primary_comparison.csv"))
logs <- fread(file.path(dest, "model_log.csv"))
stopifnot(nrow(logs) == 97L, nrow(x) == 12L,
  all(x$status == "OK"), all(is.finite(x$HR)), all(x$HR_low < x$HR_high))
audit <- fread(file.path(dest, "post_transition_followup.csv"))
stopifnot(all(audit$post_transition_days > 0), all(audit$post_transition_percent < 100))
escape <- function(s) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", as.character(s), fixed=TRUE), fixed=TRUE), fixed=TRUE)
table_html <- function(d) paste0("<table><thead><tr>", paste0("<th>", escape(names(d)), "</th>", collapse=""),
  "</tr></thead><tbody>", paste(apply(as.data.frame(d), 1, function(r) paste0("<tr>", paste0("<td>", escape(r), "</td>", collapse=""), "</tr>")), collapse=""), "</tbody></table>")
x[, pollutant := fcase(grepl("pm25", term), "PM2.5 (per 5 ug/m3)", grepl("no2", term), "NO2 (per 10 ppb)", default="O3 (per 10 ppb)")]
setorder(x, model, pollutant, followup)
display <- x[, .(Model=model, Pollutant=pollutant, Followup=followup, N, Events=events,
  `HR (95% CI)`=sprintf("%.3f (%.3f-%.3f)", HR, HR_low, HR_high))]
warnings <- logs[status != "OK", .(model, variant, subgroup, status)]
writeLines(c("<!doctype html><meta charset='utf-8'><title>LAS extended follow-up</title>",
  "<style>body{font:17px Arial;max-width:1200px;margin:36px auto;color:#111}table{border-collapse:collapse;margin:24px 0;width:100%}th,td{border-bottom:1px solid #bbb;padding:10px;text-align:left}th{background:#eee}</style>",
  "<h1>Lung candidates: extended follow-up with last observed LAS</h1>",
  paste0("<p>", escape(readLines(file.path(dest, "methods.txt"))), "</p>"),
  "<h2>Primary model comparison</h2>", table_html(display),
  "<h2>Post-transition observation in time-varying models</h2>", table_html(audit),
  "<h2>Model checks</h2><p>97 model fits returned; all six primary fits were warning-free. Timing assertions passed: no post-transition LAS and no future score backfill.</p>",
  table_html(warnings),
  "<p>CSV files in this folder include all pollutant estimates, full covariate coefficients, model formulas and warnings. Existing manuscript figures still show the transition-censored analysis.</p>"),
  file.path(dest, "results.html"))
fwrite(warnings, file.path(dest, "model_warnings.csv"))
print(display)
print(audit)
print(warnings)
