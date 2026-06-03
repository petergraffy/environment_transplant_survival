#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}
if (.Platform$OS.type == "windows") {
  r_minor <- strsplit(R.version$minor, "[.]")[[1]][[1]]
  user_lib <- Sys.getenv("R_LIBS_USER", unset = NA_character_)
  if (!is.na(user_lib) && nzchar(user_lib)) user_lib <- gsub("\\\\", "/", user_lib)
  if (is.na(user_lib) || !nzchar(user_lib)) {
    userprofile <- Sys.getenv("USERPROFILE", unset = path.expand("~"))
    user_lib <- paste0(gsub("\\\\", "/", userprofile), "/AppData/Local/R/win-library/", R.version$major, ".", r_minor)
  }
  .libPaths(c(user_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(jsonlite)
  library(readr)
  library(stringr)
  library(tibble)
})

out_dir <- file.path("data", "processed", "community")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

acs_start_year <- 2012L
acs_end_year <- 2023L
census_api_key <- Sys.getenv("CENSUS_API_KEY", unset = NA_character_)
if (is.na(census_api_key) || !nzchar(census_api_key)) {
  stop(
    "CENSUS_API_KEY is required to download ACS ZCTA covariates. ",
    "Request a Census API key and set Sys.setenv(CENSUS_API_KEY = '<key>') ",
    "or set it in the shell before running this script.",
    call. = FALSE
  )
}

acs_vars <- c(
  total_population = "B01003_001E",
  median_household_income = "B19013_001E",
  poverty_total = "B17001_001E",
  poverty_count = "B17001_002E",
  education_total_25plus = "B15003_001E",
  bachelors = "B15003_022E",
  masters = "B15003_023E",
  professional = "B15003_024E",
  doctorate = "B15003_025E",
  labor_force = "B23025_003E",
  unemployed = "B23025_005E",
  households_vehicle_total = "B08201_001E",
  households_no_vehicle = "B08201_002E",
  race_total = "B02001_001E",
  white_alone = "B02001_002E",
  black_alone = "B02001_003E",
  asian_alone = "B02001_005E",
  hispanic_total = "B03003_001E",
  hispanic = "B03003_003E"
)

safe_ratio <- function(num, den) {
  if_else(!is.na(den) & den > 0, num / den, NA_real_)
}

read_url_text <- function(url) {
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp), add = TRUE)

  method <- if (.Platform$OS.type == "windows") "wininet" else "auto"
  status <- tryCatch(
    download.file(url, tmp, quiet = TRUE, method = method),
    error = function(e) e
  )

  if (inherits(status, "error") || !file.exists(tmp) || file.info(tmp)$size == 0) {
    stop("Could not download Census API response for URL: ", sub("&key=[^&]+", "&key=<redacted>", url), call. = FALSE)
  }

  paste(readLines(tmp, warn = FALSE), collapse = "\n")
}

read_acs_zcta_year <- function(year) {
  message("Reading ACS 5-year ZCTA covariates for ", year)
  variables <- paste(c("NAME", unname(acs_vars)), collapse = ",")
  url <- paste0(
    "https://api.census.gov/data/", year, "/acs/acs5?get=",
    variables,
    "&for=zip%20code%20tabulation%20area:*",
    "&key=", census_api_key
  )

  txt <- read_url_text(url)
  if (!startsWith(trimws(txt), "[")) {
    stop("Census API did not return JSON for ", year, ". URL: ", sub("&key=[^&]+", "&key=<redacted>", url), call. = FALSE)
  }
  raw <- fromJSON(txt)
  dat <- as_tibble(raw[-1, , drop = FALSE])
  names(dat) <- make.names(raw[1, ], unique = TRUE)

  zcta_col <- grep("zip.code.tabulation.area", names(dat), value = TRUE)
  if (length(zcta_col) != 1L) {
    stop("Could not identify ZCTA column in Census API response for ", year, call. = FALSE)
  }

  names(dat)[match(unname(acs_vars), names(dat))] <- names(acs_vars)

  dat %>%
    transmute(
      zip = str_pad(.data[[zcta_col]], 5, side = "left", pad = "0"),
      community_year = year,
      total_population = as.numeric(total_population),
      median_household_income = as.numeric(median_household_income),
      pct_poverty = safe_ratio(as.numeric(poverty_count), as.numeric(poverty_total)),
      pct_bachelor_plus = safe_ratio(
        as.numeric(bachelors) + as.numeric(masters) + as.numeric(professional) + as.numeric(doctorate),
        as.numeric(education_total_25plus)
      ),
      pct_unemployed = safe_ratio(as.numeric(unemployed), as.numeric(labor_force)),
      pct_no_vehicle = safe_ratio(as.numeric(households_no_vehicle), as.numeric(households_vehicle_total)),
      pct_nonwhite = 1 - safe_ratio(as.numeric(white_alone), as.numeric(race_total)),
      pct_black = safe_ratio(as.numeric(black_alone), as.numeric(race_total)),
      pct_asian = safe_ratio(as.numeric(asian_alone), as.numeric(race_total)),
      pct_hispanic = safe_ratio(as.numeric(hispanic), as.numeric(hispanic_total))
    )
}

acs <- bind_rows(lapply(acs_start_year:acs_end_year, read_acs_zcta_year))

acs_analysis_years <- bind_rows(
  tibble(analysis_year = 2005:(acs_start_year - 1L), community_year = acs_start_year),
  tibble(analysis_year = acs_start_year:acs_end_year, community_year = acs_start_year:acs_end_year)
) %>%
  left_join(acs, by = "community_year", relationship = "many-to-many") %>%
  select(zip, analysis_year, community_year, everything())

write_csv(acs, file.path(out_dir, "zcta_acs_community_covariates_2012_2023.csv.gz"))
write_csv(acs_analysis_years, file.path(out_dir, "zcta_acs_community_covariates_2005_2023.csv.gz"))

write_csv(
  acs_analysis_years %>%
    group_by(analysis_year, community_year) %>%
    summarise(
      n_zips = n_distinct(zip),
      pct_missing_income = mean(is.na(median_household_income)),
      pct_missing_poverty = mean(is.na(pct_poverty)),
      .groups = "drop"
    ),
  file.path(out_dir, "zcta_acs_community_covariates_qc.csv")
)

message("Wrote community covariates to ", out_dir)
