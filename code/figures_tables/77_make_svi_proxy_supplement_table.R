#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(gt)
  library(readr)
  library(tibble)
})

out_dir <- file.path("output", "tables", "svi_proxy_calculation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

svi_table <- tribble(
  ~Component, ~ACS_variable_or_source, ~Calculation, ~Rank_direction, ~Use_in_proxy,
  "Poverty",
  "B17001_002E / B17001_001E",
  "Persons below poverty level divided by poverty-status denominator",
  "Higher percentage = higher vulnerability rank",
  "Percent rank within analysis year",
  "Unemployment",
  "B23025_005E / B23025_003E",
  "Unemployed civilians divided by civilian labor force",
  "Higher percentage = higher vulnerability rank",
  "Percent rank within analysis year",
  "No vehicle access",
  "B08201_002E / B08201_001E",
  "Households with no vehicle available divided by total households in vehicle table",
  "Higher percentage = higher vulnerability rank",
  "Percent rank within analysis year",
  "Non-White race",
  "1 - (B02001_002E / B02001_001E)",
  "One minus White alone population divided by total race denominator",
  "Higher percentage = higher vulnerability rank",
  "Percent rank within analysis year",
  "Low income",
  "B19013_001E",
  "Median household income",
  "Lower income = higher vulnerability rank",
  "Percent rank of negative median household income within analysis year",
  "Low educational attainment",
  "(B15003_022E + B15003_023E + B15003_024E + B15003_025E) / B15003_001E",
  "Bachelor's, master's, professional, or doctorate degree count divided by population aged 25 years or older",
  "Lower percentage with bachelor's degree or higher = higher vulnerability rank",
  "Percent rank of negative bachelor's-or-higher percentage within analysis year",
  "Overall ZCTA SVI proxy",
  "Six ranked components above",
  "Mean of poverty, unemployment, no vehicle, non-White, low-income, and low-education ranks",
  "Higher value = greater relative social vulnerability",
  "Continuous covariate added to sensitivity Cox models"
)

notes <- tribble(
  ~Note,
  "ACS variables were obtained from 5-year American Community Survey ZCTA-level tables for community years 2012 through 2023.",
  "For analysis years 2005 through 2011, 2012 ACS covariates were carried backward. For analysis years 2024 and 2025, 2023 ACS covariates were carried forward.",
  "Missing component values were median-imputed within analysis year before percentile ranks were calculated.",
  "This variable is an ACS-derived ZCTA social vulnerability proxy and is not the official CDC/ATSDR Social Vulnerability Index."
)

write_csv(svi_table, file.path(out_dir, "svi_proxy_calculation_table.csv"))
write_csv(notes, file.path(out_dir, "svi_proxy_calculation_notes.csv"))

gt_tbl <- svi_table %>%
  gt() %>%
  tab_header(title = "Supplemental Table. ACS-Derived ZCTA Social Vulnerability Proxy Calculation") %>%
  cols_label(
    Component = "Component",
    ACS_variable_or_source = "ACS variable or source",
    Calculation = "Calculation",
    Rank_direction = "Rank direction",
    Use_in_proxy = "Use in proxy"
  ) %>%
  tab_source_note("ACS indicates American Community Survey; CDC/ATSDR, Centers for Disease Control and Prevention/Agency for Toxic Substances and Disease Registry; SVI, Social Vulnerability Index; ZCTA, ZIP Code Tabulation Area.") %>%
  tab_source_note("The ZCTA SVI proxy was calculated as the mean of six within-year percentile ranks. Missing component values were median-imputed within analysis year before ranking.")

gtsave(gt_tbl, file.path(out_dir, "svi_proxy_calculation_table.html"))

message("Wrote SVI proxy supplement table to ", normalizePath(out_dir, winslash = "/"))
