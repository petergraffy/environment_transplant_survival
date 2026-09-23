#!/usr/bin/env Rscript
source("code/r_runtime.R")
ensure_user_library()
source("code/daily_pollution_inputs.R")
suppressPackageStartupMessages(library(dplyr))

scripts <- list.files("code", pattern = "^(78|79|80|81|84|86|88|90).*R$", full.names = TRUE)
invisible(lapply(scripts, parse))
summaries <- list()
for (i in seq_len(nrow(daily_pollution_releases))) {
  spec <- daily_pollution_releases[i, ]
  path <- file.path("data", "release", spec$directory)
  for (resolution in c("annual", "monthly")) {
    prior <- resolution == "annual"
    output_column <- paste0(spec$pollutant, if (prior) "_prior_" else "_interval_",
                            if (spec$pollutant == "pm25") "ug_m3" else "ppb")
    cache <- if (prior) {
      file.path("output", "prior_year_pollution_cox_svi", "cache",
                paste0(spec$pollutant, "_daily_annual_zcta.csv.gz"))
    } else {
      file.path("output", "timevarying_pollution_cox_svi", "cache",
                paste0(spec$pollutant, "_daily_derived_monthly_zcta.csv.gz"))
    }
    dat <- read_daily_pollution_aggregate(path, spec$value_column, output_column,
                                         resolution, cache)
    keys <- c("zip", "year", if (!prior) "month")
    stopifnot(setequal(unique(dat$year), 2005:2024),
              !anyDuplicated(dat[keys]), all(is.finite(dat[[output_column]])),
              all(nchar(dat$zip) == 5L),
              nrow(dat) == 33300L * 20L * if (prior) 1L else 12L)
    # Verify day weighting and leap-year inclusion against a direct daily read.
    daily <- arrow::read_parquet(tail(daily_pollution_files(path), 1L),
                                 col_select = c("zip", "date", "month", spec$value_column))
    daily <- daily[daily$zip == dat$zip[1], ]
    stopifnot(nrow(daily) == 366L, !anyDuplicated(daily$date),
              min(as.Date(daily$date)) == as.Date("2024-01-01"),
              max(as.Date(daily$date)) == daily_pollution_end_date)
    check <- dat[dat$zip == dat$zip[1] & dat$year == 2024L, ]
    expected <- if (prior) mean(daily[[spec$value_column]]) else
      vapply(check$month, function(m) mean(daily[[spec$value_column]][daily$month == m]), numeric(1))
    stopifnot(isTRUE(all.equal(as.numeric(check[[output_column]]), expected, tolerance = 1e-10)))
    summaries[[length(summaries) + 1L]] <- data.frame(
      pollutant = spec$pollutant, resolution = resolution, rows = nrow(dat),
      zctas = n_distinct(dat$zip), first_year = min(dat$year), last_year = max(dat$year),
      release = paste0("https://github.com/petergraffy/environment_transplant_survival/releases/tag/", spec$tag)
    )
    rm(dat, daily)
    gc()
  }
}
summary <- bind_rows(summaries)
dir.create("output/exposure_input_validation", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(summary, "output/exposure_input_validation/daily_exposure_cache_coverage.csv")
print(summary)
