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
  library(scales)
})

sf_use_s2(FALSE)

pm25_dir <- file.path("data", "processed", "lghap_pm25_tract_daily_cook_county_il")
o3_dir <- file.path("data", "processed", "o3_tract_daily_cook_county_il")
out_dir <- file.path("output", "figures", "cook_county_tract_pm25_o3_multiyear")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

parse_windows <- function(value) {
  if (!nzchar(value)) {
    return(list(
      `2005-2009` = 2005:2009,
      `2010-2014` = 2010:2014,
      `2015-2019` = 2015:2019,
      `2020-2021` = 2020:2021
    ))
  }
  pieces <- strsplit(value, ",")[[1]]
  out <- list()
  for (piece in trimws(pieces)) {
    bounds <- as.integer(strsplit(piece, ":", fixed = TRUE)[[1]])
    if (length(bounds) != 2L) stop("Use windows like 2005:2009,2010:2014", call. = FALSE)
    out[[sprintf("%04d-%04d", bounds[1], bounds[2])]] <- seq(bounds[1], bounds[2])
  }
  out
}

windows <- parse_windows(Sys.getenv("COOK_TRACT_MULTIYEAR_MAP_WINDOWS"))

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

window_summary <- function(years, pollutant) {
  if (pollutant == "pm25") {
    files <- as.vector(outer(years, 1:12, Vectorize(function(y, m) {
      file.path(pm25_dir, sprintf("lghap_pm25_tract_daily_%04d_%02d.csv.gz", y, m))
    })))
    value_col <- "pm25_ug_m3"
  } else {
    files <- as.vector(outer(years, 1:12, Vectorize(function(y, m) {
      file.path(o3_dir, sprintf("o3_tract_daily_%04d_%02d.csv.gz", y, m))
    })))
    value_col <- "o3_ppb"
  }

  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Missing files for ", pollutant, ":\n", paste(missing, collapse = "\n"), call. = FALSE)

  daily <- read_daily(files, value_col)
  daily[, .(
    value = mean(get(value_col), na.rm = TRUE),
    days = uniqueN(date),
    missing_values = sum(is.na(get(value_col))),
    nearest_fill_days = sum(value_source %in% c("nearest_raster_cell", "nearest_nonmissing_raster_cell"), na.rm = TRUE)
  ), by = tract_geoid]
}

hot_high_scale <- function(name, limits) {
  scale_fill_viridis_c(
    option = "magma",
    direction = 1,
    name = name,
    limits = limits,
    na.value = "grey88",
    oob = squish,
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(1.1, "in"),
      barheight = unit(5.4, "in")
    )
  )
}

plot_one <- function(summary, label, pollutant) {
  if (pollutant == "pm25") {
    title <- "PM2.5"
    legend_title <- "ug/m3"
  } else {
    title <- "Ozone"
    legend_title <- "ppb"
  }

  map_data <- merge(tracts, summary, by = "tract_geoid", all.x = TRUE)
  missing_polygons <- sum(is.na(map_data$value))
  limits <- quantile(map_data$value, c(0.02, 0.98), na.rm = TRUE)

  ggplot(map_data) +
    geom_sf(aes(fill = value), color = "white", linewidth = 0.06) +
    hot_high_scale(legend_title, limits) +
    coord_sf(datum = NA) +
    labs(
      title = sprintf("%s, %s", title, label),
      subtitle = sprintf("1,331 tracts | %d missing", missing_polygons)
    ) +
    theme_void(base_size = 24) +
    theme(
      plot.title = element_text(face = "bold", size = 36, margin = margin(b = 8)),
      plot.subtitle = element_text(size = 24, color = "grey20", margin = margin(b = 16)),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 24),
      legend.text = element_text(size = 20),
      plot.margin = margin(18, 18, 18, 18)
    )
}

all_qc <- list()
plots <- list()

for (label in names(windows)) {
  years <- windows[[label]]
  log_msg("Summarizing multi-year PM2.5 and O3 for ", label)
  pm25 <- window_summary(years, "pm25")
  o3 <- window_summary(years, "o3")

  fwrite(pm25[, .(window = label, pollutant = "pm25", tract_geoid, value, days, missing_values, nearest_fill_days)],
         file.path(out_dir, sprintf("cook_tract_pm25_multiyear_%s.csv", gsub("-", "_", label))))
  fwrite(o3[, .(window = label, pollutant = "o3", tract_geoid, value, days, missing_values, nearest_fill_days)],
         file.path(out_dir, sprintf("cook_tract_o3_multiyear_%s.csv", gsub("-", "_", label))))

  all_qc[[length(all_qc) + 1L]] <- data.table(
    window = label,
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

  panel <- plot_one(pm25, label, "pm25") + plot_one(o3, label, "o3") +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(
      title = sprintf("Cook County Tract Multi-Year Average Air Pollution, %s", label),
      subtitle = "Hot colors indicate higher concentrations",
      theme = theme(
        plot.title = element_text(face = "bold", size = 44),
        plot.subtitle = element_text(size = 30, color = "grey20")
      )
    )

  path <- file.path(out_dir, sprintf("cook_tract_pm25_o3_multiyear_%s.png", gsub("-", "_", label)))
  ggsave(path, panel, width = 18, height = 9.5, dpi = 220, bg = "white", limitsize = FALSE)
  plots[[label]] <- panel
}

qc <- rbindlist(all_qc)
fwrite(qc, file.path(out_dir, "cook_tract_pm25_o3_multiyear_map_qc.csv"))

combined <- wrap_plots(plots, ncol = 1) +
  plot_annotation(
    title = "Cook County Multi-Year PM2.5 and Ozone",
    subtitle = "Hot colors indicate higher concentrations",
    theme = theme(
      plot.title = element_text(face = "bold", size = 46),
      plot.subtitle = element_text(size = 32, color = "grey20")
    )
  )
ggsave(
  file.path(out_dir, sprintf("cook_tract_pm25_o3_multiyear_%s.png", paste(gsub("-", "_", names(windows)), collapse = "__"))),
  combined,
  width = 22,
  height = 9.5 * length(windows),
  dpi = 200,
  bg = "white",
  limitsize = FALSE
)

log_msg("Maps written to ", normalizePath(out_dir, winslash = "/"))
print(qc)
