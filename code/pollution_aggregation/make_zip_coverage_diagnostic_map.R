user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", "4.5")
if (dir.exists(user_lib)) {
  .libPaths(unique(c(user_lib, .libPaths())))
}

library(data.table)
library(ggplot2)
library(sf)

sf_use_s2(FALSE)

out_dir <- "output/zip_pollution"
map_dir <- file.path(out_dir, "maps")

zcta <- st_read("data/cache/cb_2020_us_zcta520_500k/cb_2020_us_zcta520_500k.shp", quiet = TRUE)[, c("ZCTA5CE20", "geometry")]
names(zcta)[1] <- "zip"
conus_bbox <- st_transform(
  st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326)),
  st_crs(zcta)
)
zcta <- suppressWarnings(st_crop(zcta, conus_bbox))
zcta <- st_transform(zcta, 5070)

read_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  as.data.table(read.csv(con, check.names = FALSE, colClasses = c(zip = "character")))
}

pm <- read_gz(file.path(out_dir, "pm25_2023_01.csv.gz"))[, .(zip, value_source)]
map_data <- merge(zcta, pm, by = "zip", all.x = TRUE)
map_data$coverage <- ifelse(is.na(map_data$value_source), "Missing value", "Has pollutant value")

check <- map_data |>
  st_drop_geometry() |>
  as.data.table()
check_summary <- check[, .N, by = coverage]
fwrite(check_summary, file.path(map_dir, "zip_coverage_diagnostic_summary.csv"))

p <- ggplot(map_data) +
  geom_sf(aes(fill = coverage), color = NA) +
  scale_fill_manual(
    values = c("Has pollutant value" = "#2b8cbe", "Missing value" = "#e31a1c"),
    name = NULL,
    drop = FALSE
  ) +
  coord_sf(datum = NA) +
  labs(
    title = "ZCTA Coverage Diagnostic",
    subtitle = "Blue polygons are ZIP/ZCTA areas with pollutant values; white areas are not ZCTA polygons."
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(color = "grey30", margin = margin(b = 8)),
    legend.position = "right"
  )

ggsave(
  file.path(map_dir, "zip_coverage_diagnostic.png"),
  p,
  width = 10,
  height = 6,
  dpi = 220,
  bg = "white"
)

message("Coverage diagnostic written to ", normalizePath(file.path(map_dir, "zip_coverage_diagnostic.png"), winslash = "/"))
