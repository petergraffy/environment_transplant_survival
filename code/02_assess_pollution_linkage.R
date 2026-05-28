#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(readr)
  library(sf)
  library(stringr)
  library(tibble)
})

saf_dir <- "/Users/saborpete/Library/CloudStorage/Box-Box/SAF Q2 2025"
pubsaf_dir <- file.path(saf_dir, "pubsaf2506")
supp_dir <- file.path(saf_dir, "SupplementalData2506")
pollution_dir <- "output/zip_pollution"
out_dir <- "output/linkage"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

candidate_files <- tribble(
  ~candidate_group, ~path,
  "kidney_pancreas", file.path(pubsaf_dir, "cand_kipa.sas7bdat"),
  "liver_intestine", file.path(pubsaf_dir, "cand_liin.sas7bdat"),
  "thoracic", file.path(pubsaf_dir, "cand_thor.sas7bdat")
)

clean_zip <- function(x) {
  z <- str_extract(as.character(x), "\\d{5}")
  if_else(is.na(z), NA_character_, str_pad(z, 5, side = "left", pad = "0"))
}

message("Reading candidate ZIP supplement")
candidate_zip <- read_sas(file.path(supp_dir, "canzip2506.sas7bdat")) %>%
  transmute(
    PERS_ID,
    PX_ID,
    candidate_zip = clean_zip(CAN_PERM_ZIP),
    can_wait_in_perm_zip = na_if(CAN_WAIT_IN_PERM_ZIP, "")
  ) %>%
  distinct(PERS_ID, PX_ID, .keep_all = TRUE)

message("Reading candidate SAF files")
candidates <- candidate_files %>%
  rowwise() %>%
  mutate(data = list(
    read_sas(
      path,
      col_select = any_of(c(
        "PERS_ID", "PX_ID", "WL_ORG", "CAN_LISTING_DT", "CAN_ACTIVATE_DT",
        "CAN_SOURCE", "CAN_REM_DT", "CAN_REM_CD", "CAN_DEATH_DT"
      ))
    ) %>%
      mutate(candidate_group = candidate_group)
  )) %>%
  ungroup() %>%
  select(data) %>%
  tidyr::unnest(data) %>%
  mutate(
    exposure_date = coalesce(CAN_ACTIVATE_DT, CAN_LISTING_DT),
    exposure_year = as.integer(format(exposure_date, "%Y"))
  ) %>%
  filter(exposure_year >= 2005, exposure_year <= 2024) %>%
  left_join(candidate_zip, by = c("PERS_ID", "PX_ID")) %>%
  mutate(
    has_candidate_zip = !is.na(candidate_zip) & candidate_zip != "",
    waitlist_row_id = row_number()
  )

message("Reading pollution availability")
pollution_availability <- bind_rows(
  read_csv(file.path(pollution_dir, "all_pm25_zip.csv.gz"), show_col_types = FALSE) %>%
    transmute(
      zip = clean_zip(zip),
      year,
      pollutant = "pm25",
      value = pm25_ug_m3,
      value_source
    ),
  read_csv(file.path(pollution_dir, "all_o3_zip.csv.gz"), show_col_types = FALSE) %>%
    transmute(
      zip = clean_zip(zip),
      year,
      pollutant = "o3",
      value = o3_ppb,
      value_source
    ),
  read_csv(file.path(pollution_dir, "all_no2_zip.csv.gz"), show_col_types = FALSE) %>%
    transmute(
      zip = clean_zip(zip),
      year,
      pollutant = "no2",
      value = no2,
      value_source
    )
) %>%
  filter(!is.na(zip), year >= 2005, year <= 2024) %>%
  group_by(zip, year, pollutant) %>%
  summarise(
    has_value = any(!is.na(value)),
    n_values = sum(!is.na(value)),
    has_zonal_mean = any(value_source == "zonal_mean" & !is.na(value)),
    has_nearest_fill = any(value_source == "nearest_raster_cell" & !is.na(value)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    id_cols = c(zip, year),
    names_from = pollutant,
    values_from = c(has_value, n_values, has_zonal_mean, has_nearest_fill),
    values_fill = list(
      has_value = FALSE,
      n_values = 0L,
      has_zonal_mean = FALSE,
      has_nearest_fill = FALSE
    )
  ) %>%
  mutate(
    pm25_complete_months = n_values_pm25 >= 12L,
    o3_complete_months = n_values_o3 >= 12L,
    no2_annual = n_values_no2 >= 1L,
    link_any_pollutant = has_value_pm25 | has_value_o3 | has_value_no2,
    link_all_three = pm25_complete_months & o3_complete_months & no2_annual,
    link_all_three_direct_zonal = has_zonal_mean_pm25 & has_zonal_mean_o3 & has_zonal_mean_no2,
    link_all_three_with_any_fill = link_all_three &
      (has_nearest_fill_pm25 | has_nearest_fill_o3 | has_nearest_fill_no2)
  )

linked <- candidates %>%
  left_join(
    pollution_availability,
    by = c("candidate_zip" = "zip", "exposure_year" = "year")
  ) %>%
  mutate(
    across(
      c(
        starts_with("has_value_"), starts_with("has_zonal_mean_"),
        starts_with("has_nearest_fill_"), pm25_complete_months, o3_complete_months,
        no2_annual, link_any_pollutant, link_all_three, link_all_three_direct_zonal,
        link_all_three_with_any_fill
      ),
      ~ coalesce(.x, FALSE)
    ),
    linkage_status = case_when(
      !has_candidate_zip ~ "missing_candidate_zip",
      link_all_three ~ "linked_all_three",
      link_any_pollutant ~ "partial_pollutant_linkage",
      TRUE ~ "zip_not_in_pollution_zcta"
    )
  )

registration_summary <- linked %>%
  summarise(
    unit = "waitlist_registration",
    n_total = n(),
    n_with_candidate_zip = sum(has_candidate_zip),
    n_missing_candidate_zip = sum(!has_candidate_zip),
    n_link_any_pollutant = sum(link_any_pollutant),
    n_link_all_three = sum(link_all_three),
    n_link_all_three_direct_zonal = sum(link_all_three_direct_zonal),
    n_link_all_three_with_any_fill = sum(link_all_three_with_any_fill),
    n_partial_pollutant_linkage = sum(linkage_status == "partial_pollutant_linkage"),
    n_zip_not_in_pollution_zcta = sum(linkage_status == "zip_not_in_pollution_zcta"),
    pct_link_all_three = n_link_all_three / n_total,
    pct_link_any_pollutant = n_link_any_pollutant / n_total
  )

person_summary <- linked %>%
  group_by(PERS_ID) %>%
  summarise(
    has_candidate_zip = any(has_candidate_zip),
    link_any_pollutant = any(link_any_pollutant),
    link_all_three = any(link_all_three),
    link_all_three_direct_zonal = any(link_all_three_direct_zonal),
    link_all_three_with_any_fill = any(link_all_three_with_any_fill),
    any_partial = any(linkage_status == "partial_pollutant_linkage"),
    all_zip_missing = all(!has_candidate_zip),
    all_zip_not_in_pollution_zcta = all(has_candidate_zip & !link_any_pollutant),
    .groups = "drop"
  ) %>%
  summarise(
    unit = "person",
    n_total = n(),
    n_with_candidate_zip = sum(has_candidate_zip),
    n_missing_candidate_zip = sum(all_zip_missing),
    n_link_any_pollutant = sum(link_any_pollutant),
    n_link_all_three = sum(link_all_three),
    n_link_all_three_direct_zonal = sum(link_all_three_direct_zonal),
    n_link_all_three_with_any_fill = sum(link_all_three_with_any_fill),
    n_partial_pollutant_linkage = sum(any_partial & !link_all_three),
    n_zip_not_in_pollution_zcta = sum(all_zip_not_in_pollution_zcta),
    pct_link_all_three = n_link_all_three / n_total,
    pct_link_any_pollutant = n_link_any_pollutant / n_total
  )

by_year <- linked %>%
  group_by(exposure_year) %>%
  summarise(
    n_registrations = n(),
    n_people = n_distinct(PERS_ID),
    n_missing_candidate_zip = sum(!has_candidate_zip),
    n_link_any_pollutant = sum(link_any_pollutant),
    n_link_all_three = sum(link_all_three),
    n_partial_pollutant_linkage = sum(linkage_status == "partial_pollutant_linkage"),
    n_zip_not_in_pollution_zcta = sum(linkage_status == "zip_not_in_pollution_zcta"),
    pct_link_all_three = n_link_all_three / n_registrations,
    .groups = "drop"
  )

by_organ_year <- linked %>%
  group_by(WL_ORG, exposure_year) %>%
  summarise(
    n_registrations = n(),
    n_people = n_distinct(PERS_ID),
    n_missing_candidate_zip = sum(!has_candidate_zip),
    n_link_all_three = sum(link_all_three),
    n_partial_pollutant_linkage = sum(linkage_status == "partial_pollutant_linkage"),
    n_zip_not_in_pollution_zcta = sum(linkage_status == "zip_not_in_pollution_zcta"),
    pct_link_all_three = n_link_all_three / n_registrations,
    .groups = "drop"
  )

missing_zip_examples <- linked %>%
  filter(has_candidate_zip, !link_all_three) %>%
  count(candidate_zip, exposure_year, linkage_status, sort = TRUE) %>%
  slice_head(n = 100)

zcta_path <- "data/cache/cb_2020_us_zcta520_500k/cb_2020_us_zcta520_500k.shp"
if (file.exists(zcta_path)) {
  sf_use_s2(FALSE)
  zcta <- st_read(zcta_path, quiet = TRUE) %>%
    select(candidate_zip = ZCTA5CE20)

  zcta_all <- zcta %>%
    st_drop_geometry() %>%
    transmute(
      candidate_zip = as.character(candidate_zip),
      in_2020_zcta = TRUE
    )

  conus_bbox <- st_as_sfc(st_bbox(
    c(xmin = -125, ymin = 24, xmax = -66, ymax = 50),
    crs = st_crs(zcta)
  ))
  zcta_conus <- suppressWarnings(st_crop(zcta, conus_bbox)) %>%
    st_drop_geometry() %>%
    transmute(
      candidate_zip = as.character(candidate_zip),
      in_conus_crop = TRUE
    ) %>%
    distinct()

  pollution_zips <- pollution_availability %>%
    distinct(candidate_zip = zip) %>%
    mutate(in_pollution = TRUE)

  unlinked_zip_diagnostic <- linked %>%
    filter(linkage_status == "zip_not_in_pollution_zcta") %>%
    select(PERS_ID, PX_ID, WL_ORG, exposure_year, candidate_zip) %>%
    left_join(zcta_all, by = "candidate_zip") %>%
    left_join(zcta_conus, by = "candidate_zip") %>%
    left_join(pollution_zips, by = "candidate_zip") %>%
    mutate(
      across(c(in_2020_zcta, in_conus_crop, in_pollution), ~ coalesce(.x, FALSE)),
      zip_prefix3 = substr(candidate_zip, 1, 3),
      broad_reason = case_when(
        !in_2020_zcta ~ "not_in_2020_zcta",
        in_2020_zcta & !in_conus_crop ~ "valid_2020_zcta_outside_conus_crop",
        in_conus_crop & !in_pollution ~ "in_conus_zcta_but_no_pollution",
        TRUE ~ "other"
      )
    )

  unlinked_zip_diagnostic %>%
    count(broad_reason, sort = TRUE) %>%
    mutate(pct = n / sum(n)) %>%
    write_csv(file.path(out_dir, "unlinked_zip_reason_summary.csv"))

  unlinked_zip_diagnostic %>%
    count(candidate_zip, broad_reason, sort = TRUE) %>%
    write_csv(file.path(out_dir, "unlinked_zip_diagnostic.csv"))

  unlinked_zip_diagnostic %>%
    count(broad_reason, zip_prefix3, sort = TRUE) %>%
    write_csv(file.path(out_dir, "unlinked_zip_prefix_diagnostic.csv"))
}

bind_rows(registration_summary, person_summary) %>%
  write_csv(file.path(out_dir, "pollution_linkage_overall_summary.csv"))

by_year %>%
  write_csv(file.path(out_dir, "pollution_linkage_by_year.csv"))

by_organ_year %>%
  write_csv(file.path(out_dir, "pollution_linkage_by_organ_year.csv"))

missing_zip_examples %>%
  write_csv(file.path(out_dir, "pollution_unlinked_zip_year_examples.csv"))

linked %>%
  select(
    waitlist_row_id, PERS_ID, PX_ID, WL_ORG, candidate_group, CAN_SOURCE,
    CAN_LISTING_DT, CAN_ACTIVATE_DT, exposure_year, candidate_zip,
    has_candidate_zip, link_any_pollutant, link_all_three,
    link_all_three_direct_zonal, link_all_three_with_any_fill, linkage_status
  ) %>%
  write_csv(file.path(out_dir, "candidate_pollution_linkage_flags.csv.gz"))

message("Wrote linkage outputs to ", out_dir)
