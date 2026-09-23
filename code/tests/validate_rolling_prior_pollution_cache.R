source("code/r_runtime.R")
ensure_user_library()
suppressPackageStartupMessages({ library(dplyr); library(data.table) })
source("code/rolling_prior_pollution.R")

# Independently check saved full-cohort exposures against direct daily filtering.
files <- list.files("data/cache/rolling_prior_pollution", "[.]rds$", full.names = TRUE)
for (pollutant in c("pm25", "o3")) {
  cache <- files[startsWith(basename(files), paste0(pollutant, "_"))]
  stopifnot(length(cache) > 0L)
  cache <- cache[which.max(file.info(cache)$mtime)]
  result <- readRDS(cache)
  spec <- daily_pollution_releases[daily_pollution_releases$pollutant == pollutant, ]
  value <- paste0(pollutant, "_prior_", if (pollutant == "pm25") "ug_m3" else "ppb")
  day_count <- paste0(pollutant, "_prior_days")
  stopifnot(!anyDuplicated(result[c("zip", "index_date")]),
             all(result[[day_count]] >= 0 & result[[day_count]] <= 365),
             all(is.finite(result[[value]]) == (result[[day_count]] == 365)),
             all(is.na(result[[value]][result$index_date > as.Date("2025-01-01")])))
  candidate <- result |>
    filter(index_date >= as.Date("2020-06-01"), index_date <= as.Date("2020-07-31"),
             is.finite(.data[[value]])) |>
    slice(1)
  stopifnot(nrow(candidate) == 1L)
  zip_value <- candidate$zip
  start_value <- candidate$index_date - 365
  end_value <- candidate$index_date - 1
  daily_files <- daily_pollution_files(file.path("data", "release", spec$directory))
  daily_files <- daily_files[grepl("2019|2020", daily_files)]
  actual <- arrow::open_dataset(daily_files) |>
    filter(zip == zip_value, date >= start_value, date <= end_value) |>
    select(all_of(spec$value_column)) |>
    collect()
  stopifnot(nrow(actual) == 365L,
             isTRUE(all.equal(mean(actual[[spec$value_column]]), candidate[[value]], tolerance = 1e-10)))
  cat("PASS:", pollutant, "cached coverage and direct 365-day cross-year mean including leap day.\n")
}
no2_files <- files[startsWith(basename(files), "no2_")]
no2 <- readRDS(no2_files[which.max(file.info(no2_files)$mtime)])
stopifnot(!anyDuplicated(no2[c("zip", "index_date")]),
           all(no2$no2_prior_end < no2$index_date),
           all(no2$no2_prior_months <= 12L),
           all(is.finite(no2$no2_prior_ppb) == (no2$no2_prior_months == 12L)),
           all(no2$no2_prior_annual_months[no2$index_date >= as.Date("2020-01-01")] == 0))
cat("PASS: NO2 coverage, window boundaries, and annual/monthly provenance.\n")
