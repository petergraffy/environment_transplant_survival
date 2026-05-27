user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", "4.5")
if (dir.exists(user_lib)) {
  .libPaths(unique(c(user_lib, .libPaths())))
}

library(data.table)
library(ggplot2)
library(patchwork)
library(sf)
library(viridis)

sf_use_s2(FALSE)

out_dir <- "output/zip_pollution"
map_dir <- file.path(out_dir, "maps")
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
suffix <- Sys.getenv("POLLUTION_MAP_SUFFIX", "")

zcta_path <- "data/cache/cb_2020_us_zcta520_500k/cb_2020_us_zcta520_500k.shp"
zcta <- st_read(zcta_path, quiet = TRUE)[, c("ZCTA5CE20", "geometry")]
names(zcta)[1] <- "zip"

conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326))
conus_bbox <- st_transform(conus_bbox, st_crs(zcta))
zcta <- suppressWarnings(st_crop(zcta, conus_bbox))
zcta <- st_transform(zcta, 5070)

read_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  as.data.table(read.csv(con, check.names = FALSE, colClasses = c(zip = "character")))
}

annual_mean <- function(path, value_col, map_year) {
  dt <- read_gz(path)
  dt <- dt[year == map_year]
  dt[, .(value = mean(get(value_col), na.rm = TRUE)), by = zip]
}

no2 <- read_gz(file.path(out_dir, "all_no2_zip.csv.gz"))[year == 2024, .(zip, value = no2)]
pm25 <- annual_mean(file.path(out_dir, "all_pm25_zip.csv.gz"), "pm25_ug_m3", 2023)
o3 <- annual_mean(file.path(out_dir, "all_o3_zip.csv.gz"), "o3_ppb", 2023)

datasets <- list(
  pm25 = list(data = pm25, title = "PM2.5 ZIP Mean, 2023", subtitle = "Annual mean of monthly ZIP estimates", legend = "ug/m3"),
  o3 = list(data = o3, title = "O3 ZIP Mean, 2023", subtitle = "Annual mean of monthly ZIP estimates", legend = "ppb"),
  no2 = list(data = no2, title = "NO2 ZIP Mean, 2024", subtitle = "Annual ZIP estimate", legend = "NO2")
)

plot_map <- function(data, title, subtitle, legend) {
  map_data <- merge(zcta, data, by = "zip", all.x = TRUE)
  ggplot(map_data) +
    geom_sf(aes(fill = value), color = NA) +
    scale_fill_viridis(option = "magma", direction = -1, name = legend, na.value = "grey90") +
    coord_sf(datum = NA) +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 15, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 10, color = "grey30", margin = margin(b = 8)),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      plot.margin = margin(8, 8, 8, 8)
    )
}

plots <- lapply(datasets, function(x) plot_map(x$data, x$title, x$subtitle, x$legend))

for (nm in names(plots)) {
  ggsave(
    filename = file.path(map_dir, paste0(nm, "_zip_map", suffix, ".png")),
    plot = plots[[nm]],
    width = 10,
    height = 6,
    dpi = 220,
    bg = "white"
  )
}

panel <- (plots$pm25 / plots$o3 / plots$no2) +
  plot_annotation(
    title = "ZIP-Level Air Pollution Aggregation Quick-Look Maps",
    theme = theme(plot.title = element_text(face = "bold", size = 18, margin = margin(b = 8)))
  )

ggsave(
  filename = file.path(map_dir, paste0("pollution_zip_maps_panel", suffix, ".png")),
  plot = panel,
  width = 10,
  height = 16,
  dpi = 220,
  bg = "white"
)

summary <- rbindlist(lapply(names(datasets), function(nm) {
  dt <- datasets[[nm]]$data
  data.table(
    pollutant = nm,
    mapped_zips = nrow(dt[!is.na(value)]),
    min = min(dt$value, na.rm = TRUE),
    median = median(dt$value, na.rm = TRUE),
    max = max(dt$value, na.rm = TRUE)
  )
}))
fwrite(summary, file.path(map_dir, paste0("map_value_summary", suffix, ".csv")))

map_check <- rbindlist(lapply(names(datasets), function(nm) {
  map_data <- merge(zcta, datasets[[nm]]$data, by = "zip", all.x = TRUE)
  data.table(
    pollutant = nm,
    map_polygons = nrow(map_data),
    missing_values = sum(is.na(map_data$value))
  )
}))
fwrite(map_check, file.path(map_dir, paste0("map_missingness_check", suffix, ".csv")))

message("Maps written to ", normalizePath(map_dir, winslash = "/"))
