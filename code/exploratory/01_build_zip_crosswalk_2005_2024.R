#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tibble)
  library(tidyr)
})

raw_dir <- "data/raw/hud_usps_crosswalk"
out_dir <- "data/processed/crosswalks"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

years <- 2005:2024
geographies <- c("TRACT", "COUNTY")

# HUD-USPS crosswalks are quarterly snapshots. Use Q4 for each complete year.
# For pre-2010 study years, HUD states that ZIP crosswalks are unavailable, so
# this pipeline backcasts the earliest available HUD snapshot and flags it.
hud_quarter_by_year <- tibble(
  study_year = years,
  source_year = if_else(study_year < 2010L, 2010L, study_year),
  source_quarter = if_else(study_year < 2010L, 1L, 4L),
  pre2010_backcast = study_year < 2010L
)

hud_file_stem <- function(geo, year, quarter) {
  month <- quarter * 3L
  sprintf("ZIP_%s_%02d%d", geo, month, year)
}

hud_file_candidates <- function(geo, year, quarter) {
  stem <- hud_file_stem(geo, year, quarter)
  stems <- unique(c(stem, str_to_lower(stem)))
  file.path(raw_dir, as.vector(outer(stems, c(".xlsx", ".xls", ".csv", ".txt"), paste0)))
}

hud_url_candidates <- function(geo, year, quarter) {
  stem <- hud_file_stem(geo, year, quarter)
  paste0(
    "https://www.huduser.gov/portal/datasets/usps/",
    stem,
    c(".xlsx", ".xls", ".csv", ".txt")
  )
}

find_local_hud_file <- function(geo, year, quarter) {
  candidates <- hud_file_candidates(geo, year, quarter)
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stem <- hud_file_stem(geo, year, quarter)
    file_index <- list.files(raw_dir, full.names = TRUE)
    existing <- file_index[str_to_lower(tools::file_path_sans_ext(basename(file_index))) == str_to_lower(stem)]
  }
  if (length(existing) == 0L) return(NA_character_)
  existing[[1]]
}

download_hud_file <- function(geo, year, quarter) {
  urls <- hud_url_candidates(geo, year, quarter)
  destinations <- hud_file_candidates(geo, year, quarter)

  for (i in seq_along(urls)) {
    message("Trying ", urls[[i]])
    ok <- tryCatch({
      utils::download.file(urls[[i]], destinations[[i]], mode = "wb", quiet = TRUE)
      file.exists(destinations[[i]]) && file.info(destinations[[i]])$size > 1000
    }, error = function(e) FALSE)

    if (ok) return(destinations[[i]])
    if (file.exists(destinations[[i]])) unlink(destinations[[i]])
  }

  NA_character_
}

read_hud_crosswalk <- function(path, geo) {
  ext <- tolower(tools::file_ext(path))
  dat <- if (ext %in% c("xlsx", "xls")) {
    readxl::read_excel(path, col_types = "text")
  } else {
    readr::read_delim(path, delim = ",", col_types = cols(.default = col_character()))
  }

  names(dat) <- str_to_upper(names(dat))

  required <- c("ZIP", geo, "RES_RATIO", "BUS_RATIO", "OTH_RATIO", "TOT_RATIO")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0L) {
    stop("Missing required columns in ", path, ": ", paste(missing, collapse = ", "))
  }

  dat %>%
    transmute(
      zip = str_pad(str_extract(.data$ZIP, "\\d+"), width = 5, side = "left", pad = "0"),
      geography_type = str_to_lower(geo),
      geography_id = str_pad(str_extract(.data[[geo]], "\\d+"), width = if_else(geo == "COUNTY", 5L, 11L), side = "left", pad = "0"),
      res_ratio = parse_number(.data$RES_RATIO),
      bus_ratio = parse_number(.data$BUS_RATIO),
      oth_ratio = parse_number(.data$OTH_RATIO),
      tot_ratio = parse_number(.data$TOT_RATIO),
      usps_zip_pref_city = if ("USPS_ZIP_PREF_CITY" %in% names(dat)) .data$USPS_ZIP_PREF_CITY else NA_character_,
      usps_zip_pref_state = if ("USPS_ZIP_PREF_STATE" %in% names(dat)) .data$USPS_ZIP_PREF_STATE else NA_character_
    ) %>%
    filter(!is.na(zip), !is.na(geography_id), !is.na(res_ratio))
}

load_source_crosswalk <- function(geo, source_year, source_quarter, allow_download) {
  path <- find_local_hud_file(geo, source_year, source_quarter)
  if (is.na(path) && allow_download) {
    path <- download_hud_file(geo, source_year, source_quarter)
  }
  if (is.na(path)) {
    stop(
      "No HUD-USPS file found for ZIP_", geo, " ",
      source_year, " Q", source_quarter, ".\n",
      "Place one of these files in ", raw_dir, ":\n  ",
      paste(basename(hud_file_candidates(geo, source_year, source_quarter)), collapse = "\n  ")
    )
  }

  read_hud_crosswalk(path, geo) %>%
    mutate(
      source_file = basename(path),
      source_year = source_year,
      source_quarter = source_quarter,
      .before = 1
    )
}

args <- commandArgs(trailingOnly = TRUE)
allow_download <- "--download" %in% args

manifest <- crossing(
  geography_type = geographies,
  distinct(hud_quarter_by_year, source_year, source_quarter)
) %>%
  arrange(geography_type, source_year, source_quarter)

expected_files <- manifest %>%
  mutate(
    hud_file_stem = hud_file_stem(geography_type, source_year, source_quarter),
    preferred_file = paste0(hud_file_stem, ".xlsx"),
    hud_direct_url = paste0("https://www.huduser.gov/portal/datasets/usps/", preferred_file),
    local_path = file.path(raw_dir, preferred_file)
  )

write_csv(expected_files, file.path(raw_dir, "expected_hud_usps_files_2005_2024.csv"))

all_sources <- manifest %>%
  rowwise() %>%
  mutate(data = list(load_source_crosswalk(geography_type, source_year, source_quarter, allow_download))) %>%
  ungroup() %>%
  select(-geography_type, -source_year, -source_quarter) %>%
  unnest(data)

annual_crosswalk <- hud_quarter_by_year %>%
  inner_join(
    all_sources,
    by = c("source_year", "source_quarter"),
    relationship = "many-to-many"
  ) %>%
  transmute(
    year = study_year,
    zip,
    geography_type,
    geography_id,
    res_ratio,
    bus_ratio,
    oth_ratio,
    tot_ratio,
    source = if_else(
      pre2010_backcast,
      "HUD-USPS ZIP crosswalk, 2010 Q1 backcast to pre-2010 study year",
      "HUD-USPS ZIP crosswalk"
    ),
    source_year,
    source_quarter,
    source_file,
    pre2010_backcast,
    usps_zip_pref_city,
    usps_zip_pref_state
  ) %>%
  arrange(year, geography_type, zip, desc(res_ratio), geography_id)

annual_crosswalk %>%
  filter(geography_type == "tract") %>%
  write_parquet(file.path(out_dir, "zip_to_tract_crosswalk_2005_2024.parquet"))

annual_crosswalk %>%
  filter(geography_type == "county") %>%
  write_parquet(file.path(out_dir, "zip_to_county_crosswalk_2005_2024.parquet"))

annual_crosswalk %>%
  group_by(year, geography_type, zip) %>%
  mutate(zip_res_ratio_sum = sum(res_ratio, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(year, geography_type) %>%
  summarise(
    n_rows = n(),
    n_zips = n_distinct(zip),
    n_geographies = n_distinct(geography_id),
    pct_zip_ratio_sum_near_1 = mean(abs(zip_res_ratio_sum - 1) < 0.01),
    .groups = "drop"
  ) %>%
  write_csv(file.path(out_dir, "zip_crosswalk_2005_2024_qc.csv"))

annual_crosswalk %>%
  distinct(year, source_year, source_quarter, source_file, source, pre2010_backcast) %>%
  arrange(year, source_file) %>%
  write_csv(file.path(out_dir, "zip_crosswalk_2005_2024_sources.csv"))

message("Wrote crosswalk outputs to ", out_dir)
