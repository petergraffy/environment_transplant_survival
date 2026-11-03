source("code/r_runtime.R")
ensure_user_library()
library(data.table)
source("code/daily_pollution_inputs.R")
cohort <- data.table(PX_ID = 1:4, PERS_ID = 1:4, waitlist_row_id = 1:4,
                     WL_ORG = c("HR", "KI", "LI", "LU"), candidate_zip = "01234",
                     index_date = as.Date(c("2020-01-31", "2020-02-28", "2023-12-31", "2025-01-01")),
                     observed_end_date = as.Date(c("2020-03-02", "2020-04-12", "2025-01-04", "2025-02-01")),
                     adverse_event = c(1L, 0L, 1L, 0L), age = 50, sex = "F", race = "WHITE",
                     zcta_svi_proxy = .5, listing_center = "A", listing_year = factor(c(2020,2020,2023,2025)),
                     organ_score = 3, kidney_dialysis_years = 2, kidney_no_dialysis_time = 0L,
                     kidney_diabetes = 1L, age_group = "40-64", race_group = "White",
                     multi_organ_candidate = "Single-organ candidate", unused_lab = 99)
extract <- function(path) {
  statements <- parse(path)
  item <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                   identical(x[[2]], as.name("make_month_intervals")), as.list(statements))
  stopifnot(length(item) == 1L)
  env <- new.env(parent = globalenv())
  eval(item[[1]], env)
  env$make_month_intervals
}
paths <- list.files("code", pattern = "^(79|84|86|88).*R$", full.names = TRUE)
for (path in paths) {
  previous <- extract(file.path("output/revision_runs/rolling_20260922/code_at_start", basename(path)))
  current <- extract(path)
  old <- previous(copy(cohort), as.Date("2024-12-31"))
  new <- current(copy(cohort), as.Date("2024-12-31"))
  stopifnot(!"unused_lab" %in% names(new), identical(old[, names(new), with = FALSE], new))
}
cat("PASS: monthly intervals, outcomes, and all retained model variables are identical after trimming unused columns.\n")
