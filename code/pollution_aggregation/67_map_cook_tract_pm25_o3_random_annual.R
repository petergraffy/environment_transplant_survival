#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}
local_appdata <- Sys.getenv("LOCALAPPDATA", unset = NA_character_)
userprofile <- Sys.getenv("USERPROFILE", unset = NA_character_)
r_minor <- strsplit(R.version$minor, "[.]", fixed = FALSE)[[1]][[1]]
candidate_libs <- character()
if (!is.na(local_appdata) && nzchar(local_appdata)) {
  candidate_libs <- c(candidate_libs, file.path(local_appdata, "R", "win-library", paste0(R.version$major, ".", r_minor)))
}
if (!is.na(userprofile) && nzchar(userprofile)) {
  candidate_libs <- c(candidate_libs, file.path(userprofile, "AppData", "Local", "R", "win-library", paste0(R.version$major, ".", r_minor)))
}
candidate_libs <- candidate_libs[dir.exists(candidate_libs)]
if (length(candidate_libs)) .libPaths(c(candidate_libs, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(sf)
  library(viridis)
})

sf_use_s2(FALSE)

pm25_dir <- file.path("data", "processed", "lghap_pm25_tract_daily_cook_county_il")
o3_dir <- file.path("data", "processed", "o3_tract_daily_cook_county_il")
out_dir <- file.path("output", "figures", "cook_county_tract_pm25_o3_random_annual")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

years_env <- Sys.getenv("COOK_TRACT_ANNUAL_MAP_YEARS")
if (nzchar(years_env)) {
  target_years <- as.integer(strsplit(years_env, ",")[[1]])
} else {
  set.seed(as.integer(Sys.getenv("COOK_TRACT_ANNUAL_MAP_SEED", "20260624")))
  target_years <- sort(sample(2005:2021, 4))
}

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

tract_path <- file.path("data", "cache", "cb_2020_us_tract_500k", "cb_2020_us_tract_500k.shp")
if (!file.exists(tract_path)) stop("Missing cached 2020 tract shapefile: ", tract_path, call. = FALSE)

log_msg("Reading 2020 Cook County tract boundaries")
tracts <- st_read(tract_path, quiet = TRUE)
tracts <- tracts[tracts$STATEFP == "17" & tracts$COUNTYFP == "031", c("GEOID", "geometry")]
names(tracts)[1] <- "tract_geoid"
tracts$tract_geoid <- as.character(tracts$tract_geoid)
tracts <- st_transform(st_make_valid(tracts), 26916)

read_daily <- function(files, value_col) {
  rbindlist(lapply(files, function(path) {
    log_msg("Reading ", basename(path))
    dt <- as.data.table(read.csv(
      gzfile(path),
      colClasses = c(
        tract_geoid = "character",
        date = "character",
        pollutant = "character",
        year = "integer",
        month = "integer",
        temporal_resolution = "character",
        pm25_ug_m3 = "numeric",
        o3_ppb = "numeric",
        value_source = "character",
        fill_distance_m = "numeric",
        fill_cell = "numeric"
      )
    ))
    dt <- dt[, c("tract_geoid", "date", value_col, "value_source"), with = FALSE]
    dt[, `:=`(tract_geoid = sprintf("%011s", tract_geoid), date = as.Date(date))]
    dt
  }), fill = TRUE)
}

annual_summary <- function(year, pollutant) {
  if (pollutant == "pm25") {
    files <- file.path(pm25_dir, sprintf("lghap_pm25_tract_daily_%04d_%02d.csv.gz", year, 1:12))
    value_col <- "pm25_ug_m3"
  } else {
    files <- file.path(o3_dir, sprintf("o3_tract_daily_%04d_%02d.csv.gz", year, 1:12))
    value_col <- "o3_ppb"
  }

  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Missing files for ", pollutant, " ", year, ":\n", paste(missing, collapse = "\n"), call. = FALSE)

  daily <- read_daily(files, value_col)
  daily[, .(
    value = mean(get(value_col), na.rm = TRUE),
    days = uniqueN(date),
    missing_values = sum(is.na(get(value_col))),
    nearest_fill_days = sum(value_source %in% c("nearest_raster_cell", "nearest_nonmissing_raster_cell"), na.rm = TRUE)
  ), by = tract_geoid]
}

plot_one <- function(summary, year, pollutant) {
  if (pollutant == "pm25") {
    title <- "PM2.5"
    legend_title <- "ug/m3"
    option <- "magma"
  } else {
    title <- "Ozone"
    legend_title <- "ppb"
    option <- "plasma"
  }

  map_data <- merge(tracts, summary, by = "tract_geoid", all.x = TRUE)
  missing_polygons <- sum(is.na(map_data$value))
  limits <- quantile(map_data$value, c(0.02, 0.98), na.rm = TRUE)

  ggplot(map_data) +
    geom_sf(aes(fill = pmin(pmax(value, limits[1]), limits[2])), color = "white", linewidth = 0.04) +
    scale_fill_viridis(
      option = option,
      direction = -1,
      name = legend_title,
      limits = limits,
      na.value = "grey88"
    ) +
    coord_sf(datum = NA) +
    labs(
      title = sprintf("%s, %s", title, year),
      subtitle = sprintf("%d tracts, %d missing", nrow(map_data), missing_polygons)
    ) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 11, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 8, color = "grey35", margin = margin(b = 5)),
      legend.position = "right",
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      plot.margin = margin(5, 5, 5, 5)
    )
}

all_qc <- list()
plots <- list()

for (year in target_years) {
  log_msg("Summarizing annual PM2.5 and O3 for ", year)
  pm25 <- annual_summary(year, "pm25")
  o3 <- annual_summary(year, "o3")

  fwrite(pm25[, .(year = year, pollutant = "pm25", tract_geoid, value, days, missing_values, nearest_fill_days)],
         file.path(out_dir, sprintf("cook_tract_pm25_annual_%04d.csv", year)))
  fwrite(o3[, .(year = year, pollutant = "o3", tract_geoid, value, days, missing_values, nearest_fill_days)],
         file.path(out_dir, sprintf("cook_tract_o3_annual_%04d.csv", year)))

  all_qc[[length(all_qc) + 1L]] <- data.table(
    year = year,
    pollutant = c("pm25", "o3"),
    tracts = c(nrow(pm25), nrow(o3)),
    min_days = c(min(pm25$days), min(o3$days)),
    max_days = c(max(pm25$days), max(o3$days)),
    missing_values = c(sum(pm25$missing_values), sum(o3$missing_values)),
    missing_tracts = c(sum(is.na(pm25$value)), sum(is.na(o3$value))),
    mean_value = c(mean(pm25$value, na.rm = TRUE), mean(o3$value, na.rm = TRUE)),
    min_value = c(min(pm25$value, na.rm = TRUE), min(o3$value, na.rm = TRUE)),
    max_value = c(max(pm25$value, na.rm = TRUE), max(o3$value, na.rm = TRUE))
  )

  year_panel <- plot_one(pm25, year, "pm25") + plot_one(o3, year, "o3") +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(
      title = sprintf("Cook County Census Tract Annual Average Air Pollution, %s", year),
      subtitle = "Daily tract estimates averaged over the calendar year",
      theme = theme(plot.title = element_text(face = "bold", size = 15))
    )

  year_path <- file.path(out_dir, sprintf("cook_tract_pm25_o3_annual_%04d.png", year))
  ggsave(year_path, year_panel, width = 12, height = 6, dpi = 240, bg = "white")
  plots[[as.character(year)]] <- year_panel
}

qc <- rbindlist(all_qc)
fwrite(qc, file.path(out_dir, "cook_tract_pm25_o3_random_annual_map_qc.csv"))

combined <- wrap_plots(plots, ncol = 1) +
  plot_annotation(
    title = sprintf("Random Annual Cook County Tract PM2.5 and Ozone Maps: %s", paste(target_years, collapse = ", ")),
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )
ggsave(
  file.path(out_dir, sprintf("cook_tract_pm25_o3_random_annual_%s.png", paste(target_years, collapse = "_"))),
  combined,
  width = 12,
  height = 6 * length(target_years),
  dpi = 220,
  bg = "white",
  limitsize = FALSE
)

log_msg("Maps written to ", normalizePath(out_dir, winslash = "/"))
print(qc)
