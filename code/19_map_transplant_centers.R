#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(haven)
  library(readr)
  library(sf)
  library(stringr)
  library(tibble)
})

saf_paths_source <- file.path("code", "saf_paths.R")
if (!file.exists(saf_paths_source)) saf_paths_source <- "saf_paths.R"
source(saf_paths_source)
saf_paths <- get_saf_paths()
assert_saf_files(saf_paths, include_stathist = FALSE)

pubsaf_dir <- saf_paths$pubsaf_dir
out_dir <- file.path("output", "center_maps")
fig_dir <- file.path("output", "figures", "center_maps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

analysis_start_year <- 2006L
analysis_end_year <- 2023L
target_organs <- c("HR", "KI", "LI", "LU")

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

clean_zip <- function(x) {
  z <- str_extract(as.character(x), "[0-9]+")
  z <- if_else(nchar(z) > 5L, substr(z, 1, 5), z)
  if_else(is.na(z), NA_character_, str_pad(z, 5, side = "left", pad = "0"))
}

organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")

candidate_files <- tribble(
  ~candidate_group, ~path,
  "kidney_pancreas", file.path(pubsaf_dir, "cand_kipa.sas7bdat"),
  "liver_intestine", file.path(pubsaf_dir, "cand_liin.sas7bdat"),
  "thoracic", file.path(pubsaf_dir, "cand_thor.sas7bdat")
)

candidate_cols <- c(
  "PERS_ID", "PX_ID", "WL_ORG", "CAN_LISTING_DT", "CAN_ACTIVATE_DT",
  "CAN_SOURCE", "CAN_REM_DT", "CAN_DEATH_DT", "CAN_ENDWLFU",
  "CAN_LISTING_CTR_CD"
)

log_msg("Reading candidate listing centers")
center_counts <- candidate_files %>%
  rowwise() %>%
  mutate(data = list(read_sas(path, col_select = any_of(candidate_cols)) %>% mutate(candidate_group = candidate_group))) %>%
  ungroup() %>%
  select(data) %>%
  tidyr::unnest(data) %>%
  filter(WL_ORG %in% target_organs) %>%
  mutate(
    index_date = coalesce(CAN_ACTIVATE_DT, CAN_LISTING_DT),
    index_year = as.integer(format(index_date, "%Y")),
    event_date = coalesce(CAN_REM_DT, CAN_DEATH_DT, CAN_ENDWLFU),
    center_code = as.character(CAN_LISTING_CTR_CD)
  ) %>%
  filter(
    index_year >= analysis_start_year,
    index_year <= analysis_end_year,
    CAN_SOURCE %in% c("A", "R"),
    !is.na(index_date),
    !is.na(event_date),
    !is.na(center_code),
    center_code != ""
  ) %>%
  count(center_code, organ = WL_ORG, name = "n_candidates") %>%
  mutate(organ_label = recode(organ, !!!organ_labels))

center_totals <- center_counts %>%
  group_by(center_code) %>%
  summarise(
    total_candidates = sum(n_candidates),
    dominant_organ = organ_label[which.max(n_candidates)],
    n_organs = n_distinct(organ),
    .groups = "drop"
  )

log_msg("Reading SAF institution table")
institution <- read_sas(file.path(pubsaf_dir, "institution.sas7bdat")) %>%
  transmute(
    center_code = as.character(CTR_CD),
    center_type = CTR_TY,
    optn_region = REGION,
    center_name = ENTIRE_NAME,
    city = PRIMARY_CITY,
    state = PRIMARY_STATE,
    zip = clean_zip(PRIMARY_ZIP),
    country = PRIMARY_CTRY
  ) %>%
  filter(center_type == "TX1", country %in% c("US", "USA", "United States", "", NA_character_)) %>%
  distinct(center_code, .keep_all = TRUE)

log_msg("Reading ZCTA centroids")
zcta_shp <- file.path("data", "cache", "cb_2020_us_zcta520_500k", "cb_2020_us_zcta520_500k.shp")
if (!file.exists(zcta_shp)) stop("Missing cached ZCTA shapefile: ", zcta_shp, call. = FALSE)
zcta_centroids <- st_read(zcta_shp, quiet = TRUE) %>%
  st_transform(4326) %>%
  mutate(zip = ZCTA5CE20) %>%
  st_point_on_surface() %>%
  mutate(
    lon = st_coordinates(.)[, 1],
    lat = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  select(zip, lon, lat)

city_layer <- NULL
if (requireNamespace("maps", quietly = TRUE)) {
  city_layer <- maps::us.cities %>%
    as_tibble() %>%
    transmute(
      city_key = str_to_lower(str_squish(str_remove(name, paste0(" ", country.etc, "$")))),
      state = country.etc,
      city_lon = long,
      city_lat = lat,
      pop
    ) %>%
    group_by(city_key, state) %>%
    slice_max(pop, n = 1, with_ties = FALSE) %>%
    ungroup()
}

add_location <- function(data) {
  data <- data %>%
    mutate(city_key = str_to_lower(str_squish(city)))
  if (!is.null(city_layer)) {
    data <- data %>%
      left_join(city_layer, by = c("city_key", "state")) %>%
      mutate(
        location_source = case_when(
          !is.na(lon) & !is.na(lat) ~ "zip_zcta",
          is.na(lon) & is.na(lat) & !is.na(city_lon) & !is.na(city_lat) ~ "city_state",
          TRUE ~ NA_character_
        ),
        lon = coalesce(lon, city_lon),
        lat = coalesce(lat, city_lat)
      ) %>%
      select(-city_lon, -city_lat, -pop)
  } else {
    data <- data %>% mutate(location_source = if_else(!is.na(lon) & !is.na(lat), "zip_zcta", NA_character_))
  }
  data %>% select(-city_key)
}

centers_all <- center_totals %>%
  left_join(institution, by = "center_code") %>%
  left_join(zcta_centroids, by = "zip") %>%
  add_location()

center_organs_all <- center_counts %>%
  left_join(institution, by = "center_code") %>%
  left_join(zcta_centroids, by = "zip") %>%
  add_location()

centers <- centers_all %>%
  filter(!is.na(lon), !is.na(lat), lon >= -125, lon <= -66, lat >= 24, lat <= 50)

center_organs <- center_organs_all %>%
  filter(!is.na(lon), !is.na(lat), lon >= -125, lon <= -66, lat >= 24, lat <= 50)

unmapped <- centers_all %>%
  filter(is.na(lon) | is.na(lat) | lon < -125 | lon > -66 | lat < 24 | lat > 50)

write_csv(centers, file.path(out_dir, "transplant_centers_mapped.csv"))
write_csv(center_organs, file.path(out_dir, "transplant_centers_by_organ_mapped.csv"))
write_csv(unmapped, file.path(out_dir, "transplant_centers_unmapped.csv"))

state_layer <- NULL
if (requireNamespace("maps", quietly = TRUE)) {
  state_layer <- ggplot2::map_data("state")
}

base_map <- function() {
  p <- ggplot()
  if (!is.null(state_layer)) {
    p <- p + geom_polygon(
      data = state_layer,
      aes(x = long, y = lat, group = group),
      fill = "#f7f7f4",
      color = "#cfcfca",
      linewidth = 0.25
    )
  }
  p +
    coord_quickmap(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
    theme_void(base_size = 11) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      plot.title = element_text(face = "bold", margin = margin(b = 4)),
      plot.subtitle = element_text(margin = margin(b = 8)),
      strip.text = element_text(face = "bold")
    )
}

p_all <- base_map() +
  geom_point(
    data = centers,
    aes(x = lon, y = lat, size = total_candidates, color = dominant_organ),
    alpha = 0.72
  ) +
  scale_size_continuous(range = c(1.4, 8), labels = scales::label_comma()) +
  scale_color_manual(values = c(Heart = "#b6403a", Kidney = "#1f6f8b", Liver = "#5f8f3a", Lung = "#7b5ab6")) +
  labs(
    title = "US transplant centers represented in the analytic cohort",
    subtitle = paste0(analysis_start_year, "-", analysis_end_year, " active waitlist candidates; point size is candidate volume, color is dominant organ"),
    size = "Candidates",
    color = "Dominant organ"
  )

p_organ <- base_map() +
  geom_point(
    data = center_organs,
    aes(x = lon, y = lat, size = n_candidates),
    color = "#1f6f8b",
    alpha = 0.72
  ) +
  facet_wrap(~organ_label, ncol = 2) +
  scale_size_continuous(range = c(1.1, 6.5), labels = scales::label_comma()) +
  labs(
    title = "Transplant center locations by organ",
    subtitle = paste0(analysis_start_year, "-", analysis_end_year, " active waitlist candidates; point size is organ-specific candidate volume"),
    size = "Candidates"
  )

ggsave(file.path(fig_dir, "transplant_centers_all_organs_map.png"), p_all, width = 11, height = 6.8, dpi = 320)
ggsave(file.path(fig_dir, "transplant_centers_by_organ_map.png"), p_organ, width = 11, height = 8.2, dpi = 320)

message("Mapped ", nrow(centers), " centers; wrote maps to ", fig_dir)
