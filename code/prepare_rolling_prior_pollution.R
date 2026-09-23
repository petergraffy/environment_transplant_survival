#!/usr/bin/env Rscript
source("code/r_runtime.R")
ensure_user_library()
suppressPackageStartupMessages({ library(data.table); library(dplyr) })
source("code/rolling_prior_pollution.R")
path <- "output/primary_waitlist_period_pollution_cox/primary_waitlist_period_pollution_analysis_dataset.csv.gz"
cohort <- readr::read_csv(path, col_select = c(WL_ORG, candidate_zip, index_date),
                          col_types = readr::cols(candidate_zip = readr::col_character()), show_col_types = FALSE) |>
  dplyr::select(WL_ORG, candidate_zip, index_date) |>
  dplyr::mutate(index_date = as.Date(index_date))
exposures <- read_rolling_prior_pollution(cohort, annual_fallback = TRUE)
out <- "output/exposure_input_validation"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
summary <- dplyr::bind_rows(lapply(names(exposures), function(p) {
  values <- exposures[[p]]
  term <- paste0(p, "_prior_", if (p == "pm25") "ug_m3" else "ppb")
  joined <- dplyr::left_join(cohort, values, by = c("candidate_zip" = "zip", "index_date"))
  joined |>
    dplyr::mutate(listing_year = as.integer(format(index_date, "%Y")),
                   complete = is.finite(.data[[term]])) |>
    dplyr::group_by(WL_ORG, listing_year) |>
    dplyr::summarise(candidates = dplyr::n(), complete_windows = sum(complete),
                     incomplete_windows = sum(!complete), .groups = "drop") |>
    dplyr::mutate(pollutant = p)
}))
readr::write_csv(summary, file.path(out, "rolling_prior_exposure_coverage_by_organ_year.csv"))
no2 <- exposures$no2
readr::write_csv(as.data.frame(table(no2$no2_prior_annual_months, useNA = "ifany")),
                  file.path(out, "rolling_no2_annual_approximation_months.csv"))
print(summary |> group_by(pollutant) |> summarise(candidates = sum(candidates), complete_windows = sum(complete_windows)))
