#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

old_root <- "/Users/saborpete/Library/CloudStorage/Box-Box/SAF Q2 2025"
new_root <- "/Users/saborpete/Library/CloudStorage/Box-Box/Q1 2026 SAF"
out_dir <- "output/saf_dictionary_comparison"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

old_pubsaf <- file.path(old_root, "pubsaf2506")
new_pubsaf <- file.path(new_root, "pubsaf2603")

old_supp_dirs <- c(
  file.path(old_root, "SupplementalData2506"),
  file.path(old_root, "Supplement25061"),
  file.path(old_root, "Supplement25061", "ThoracicRegistration")
)
new_supp_dirs <- c(
  file.path(new_root, "Supplement2603"),
  file.path(new_root, "Supplement2603", "Thoracic Registration")
)

clean_table_name <- function(path) {
  nm <- tools::file_path_sans_ext(basename(path))
  nm <- str_replace(nm, "2506$", "")
  nm <- str_replace(nm, "2603$", "")
  nm
}

list_sas_files <- function(dirs, release, section) {
  tibble(path = unlist(lapply(dirs, function(d) {
    if (!dir.exists(d)) return(character())
    list.files(d, pattern = "\\.sas7bdat$", full.names = TRUE, recursive = FALSE)
  }))) %>%
    mutate(
      release = release,
      section = section,
      file_name = basename(path),
      table = clean_table_name(path),
      rel_path = str_remove(path, paste0("^", fixed(if_else(release == "old_2506", old_root, new_root)), "/?"))
    ) %>%
    select(release, section, table, file_name, rel_path, path) %>%
    arrange(section, table)
}

extract_metadata <- function(path) {
  dat <- read_sas(path, n_max = 0)
  tibble(
    variable_order = seq_along(dat),
    variable = names(dat),
    storage_class = map_chr(dat, ~ paste(class(.x), collapse = "|")),
    label = map_chr(dat, ~ {
      lab <- attr(.x, "label", exact = TRUE)
      if (is.null(lab)) "" else as.character(lab)
    }),
    sas_format = map_chr(dat, ~ {
      fmt <- attr(.x, "format.sas", exact = TRUE)
      if (is.null(fmt)) "" else as.character(fmt)
    })
  )
}

safe_metadata <- function(path) {
  tryCatch(
    extract_metadata(path),
    error = function(e) {
      tibble(
        variable_order = NA_integer_,
        variable = NA_character_,
        storage_class = NA_character_,
        label = paste0("READ_ERROR: ", conditionMessage(e)),
        sas_format = NA_character_
      )
    }
  )
}

message("Listing SAF SAS files")
inventory <- bind_rows(
  list_sas_files(old_pubsaf, "old_2506", "pubsaf"),
  list_sas_files(new_pubsaf, "new_2603", "pubsaf"),
  list_sas_files(old_supp_dirs, "old_2506", "supplement"),
  list_sas_files(new_supp_dirs, "new_2603", "supplement")
)
write_csv(inventory %>% select(-path), file.path(out_dir, "saf_file_inventory.csv"))

inventory_comparison <- inventory %>%
  count(section, table, release, name = "n_files") %>%
  pivot_wider(names_from = release, values_from = n_files, values_fill = 0) %>%
  mutate(
    table_status = case_when(
      old_2506 > 0 & new_2603 > 0 ~ "present_in_both",
      old_2506 > 0 & new_2603 == 0 ~ "old_only",
      old_2506 == 0 & new_2603 > 0 ~ "new_only",
      TRUE ~ "unknown"
    )
  ) %>%
  arrange(section, table_status, table)
write_csv(inventory_comparison, file.path(out_dir, "saf_table_inventory_comparison.csv"))

message("Reading zero-row SAS metadata")
metadata <- inventory %>%
  mutate(metadata = map(path, safe_metadata)) %>%
  select(release, section, table, file_name, rel_path, metadata) %>%
  unnest(metadata)
write_csv(metadata, file.path(out_dir, "saf_variable_metadata_long.csv"))

variable_presence <- metadata %>%
  filter(!is.na(variable)) %>%
  distinct(section, table, variable, release) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = release, values_from = present, values_fill = FALSE) %>%
  mutate(
    variable_status = case_when(
      old_2506 & new_2603 ~ "present_in_both",
      old_2506 & !new_2603 ~ "old_only",
      !old_2506 & new_2603 ~ "new_only",
      TRUE ~ "unknown"
    )
  ) %>%
  arrange(section, table, variable_status, variable)
write_csv(variable_presence, file.path(out_dir, "saf_variable_presence_comparison.csv"))

variable_changes <- metadata %>%
  filter(!is.na(variable)) %>%
  select(release, section, table, variable, variable_order, storage_class, label, sas_format) %>%
  pivot_wider(
    names_from = release,
    values_from = c(variable_order, storage_class, label, sas_format),
    names_sep = "__"
  ) %>%
  filter(!is.na(variable_order__old_2506), !is.na(variable_order__new_2603)) %>%
  mutate(
    order_changed = variable_order__old_2506 != variable_order__new_2603,
    storage_class_changed = coalesce(storage_class__old_2506, "") != coalesce(storage_class__new_2603, ""),
    label_changed = coalesce(label__old_2506, "") != coalesce(label__new_2603, ""),
    sas_format_changed = coalesce(sas_format__old_2506, "") != coalesce(sas_format__new_2603, "")
  ) %>%
  filter(order_changed | storage_class_changed | label_changed | sas_format_changed) %>%
  arrange(section, table, variable)
write_csv(variable_changes, file.path(out_dir, "saf_variable_metadata_changes.csv"))

table_summary <- metadata %>%
  filter(!is.na(variable)) %>%
  count(release, section, table, name = "n_variables") %>%
  pivot_wider(names_from = release, values_from = n_variables) %>%
  mutate(
    n_variables_old_2506 = coalesce(old_2506, 0L),
    n_variables_new_2603 = coalesce(new_2603, 0L),
    n_variable_delta = n_variables_new_2603 - n_variables_old_2506
  ) %>%
  select(section, table, n_variables_old_2506, n_variables_new_2603, n_variable_delta) %>%
  arrange(section, desc(abs(n_variable_delta)), table)
write_csv(table_summary, file.path(out_dir, "saf_table_variable_count_comparison.csv"))

analysis_vars <- tribble(
  ~table, ~variable,
  "cand_kipa", "PERS_ID",
  "cand_kipa", "PX_ID",
  "cand_kipa", "WL_ORG",
  "cand_kipa", "CAN_LISTING_DT",
  "cand_kipa", "CAN_ACTIVATE_DT",
  "cand_kipa", "CAN_SOURCE",
  "cand_kipa", "CAN_REM_DT",
  "cand_kipa", "CAN_REM_CD",
  "cand_kipa", "CAN_DEATH_DT",
  "cand_kipa", "CAN_ENDWLFU",
  "cand_kipa", "CAN_AGE_AT_LISTING",
  "cand_kipa", "CAN_GENDER",
  "cand_kipa", "CAN_RACE_SRTR",
  "cand_kipa", "CAN_ETHNICITY_SRTR",
  "cand_kipa", "CAN_BMI",
  "cand_liin", "PERS_ID",
  "cand_liin", "PX_ID",
  "cand_liin", "WL_ORG",
  "cand_liin", "CAN_LISTING_DT",
  "cand_liin", "CAN_ACTIVATE_DT",
  "cand_liin", "CAN_SOURCE",
  "cand_liin", "CAN_REM_DT",
  "cand_liin", "CAN_REM_CD",
  "cand_liin", "CAN_DEATH_DT",
  "cand_liin", "CAN_ENDWLFU",
  "cand_liin", "CAN_AGE_AT_LISTING",
  "cand_liin", "CAN_GENDER",
  "cand_liin", "CAN_RACE_SRTR",
  "cand_liin", "CAN_ETHNICITY_SRTR",
  "cand_liin", "CAN_BMI",
  "cand_thor", "PERS_ID",
  "cand_thor", "PX_ID",
  "cand_thor", "WL_ORG",
  "cand_thor", "CAN_LISTING_DT",
  "cand_thor", "CAN_ACTIVATE_DT",
  "cand_thor", "CAN_SOURCE",
  "cand_thor", "CAN_REM_DT",
  "cand_thor", "CAN_REM_CD",
  "cand_thor", "CAN_DEATH_DT",
  "cand_thor", "CAN_ENDWLFU",
  "cand_thor", "CAN_AGE_AT_LISTING",
  "cand_thor", "CAN_GENDER",
  "cand_thor", "CAN_RACE_SRTR",
  "cand_thor", "CAN_ETHNICITY_SRTR",
  "cand_thor", "CAN_BMI",
  "canzip", "PERS_ID",
  "canzip", "PX_ID",
  "canzip", "CAN_PERM_ZIP",
  "stathist_kipa", "PX_ID",
  "stathist_kipa", "WL_ORG",
  "stathist_kipa", "CANHX_BEGIN_DT",
  "stathist_kipa", "CANHX_END_DT",
  "stathist_kipa", "CANHX_STAT_CD",
  "stathist_liin", "PX_ID",
  "stathist_liin", "WL_ORG",
  "stathist_liin", "CANHX_BEGIN_DT",
  "stathist_liin", "CANHX_END_DT",
  "stathist_liin", "CANHX_STAT_CD",
  "stathist_thor", "PX_ID",
  "stathist_thor", "WL_ORG",
  "stathist_thor", "CANHX_BEGIN_DT",
  "stathist_thor", "CANHX_END_DT",
  "stathist_thor", "CANHX_STAT_CD"
)

analysis_variable_check <- analysis_vars %>%
  left_join(
    metadata %>%
      filter(!is.na(variable)) %>%
      select(release, section, table, variable, storage_class, label, sas_format),
    by = c("table", "variable")
  ) %>%
  mutate(present = !is.na(release)) %>%
  select(table, variable, release, section, present, storage_class, label, sas_format) %>%
  arrange(table, variable, release)
write_csv(analysis_variable_check, file.path(out_dir, "analysis_variable_check.csv"))

message("Wrote SAF dictionary comparison outputs to ", out_dir)
