user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", "4.5")
if (dir.exists(user_lib)) {
  .libPaths(unique(c(user_lib, .libPaths())))
}

library(data.table)
library(sf)
library(terra)
library(viridis)

raw_map_dir <- "output/zip_pollution/maps/raw_grids"
dir.create(raw_map_dir, recursive = TRUE, showWarnings = FALSE)

raw_files <- data.table(
  pollutant = c("pm25", "o3", "no2"),
  label = c("PM2.5 raw grid, January 2023", "O3 raw grid, January 2023", "NO2 raw grid, annual 2024"),
  file = c(
    "C:/Users/Peter Graffy/Box/PM2.5/Monthly/Monthly/2023/V6GL02.04.CNNPM25.NA.202301-202301.nc",
    "C:/Users/Peter Graffy/Box/O3/CONUS_O3_MDA8_p01_monthly_mean_2023.nc",
    "C:/Users/Peter Graffy/Box/NO2/Anenberg et al 2022 Global NO2/annual_mean_tropomi_lur_conus_surface_no2_2024.v1.02.nc"
  ),
  layer = c(1L, 1L, 1L),
  units = c("ug/m3", "ppb", "NO2")
)

conus_ext <- ext(-125, -66, 24, 50)
zcta <- st_read("data/cache/cb_2020_us_zcta520_500k/cb_2020_us_zcta520_500k.shp", quiet = TRUE)
conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326))
conus_bbox <- st_transform(conus_bbox, st_crs(zcta))
zcta <- suppressWarnings(st_crop(zcta, conus_bbox))
zcta_vect <- vect(zcta)
summary_rows <- vector("list", nrow(raw_files))

for (i in seq_len(nrow(raw_files))) {
  meta <- raw_files[i]
  r <- rast(meta$file)[[meta$layer]]
  if (is.na(crs(r))) crs(r) <- "EPSG:4326"
  r <- crop(r, conus_ext)
  if (meta$pollutant == "pm25") {
    r[r <= -900] <- NA
  }

  zmask <- rasterize(project(zcta_vect, crs(r)), r, field = 1, touches = TRUE)
  na_count <- as.numeric(global(is.na(r), "sum", na.rm = FALSE)[1, 1])
  cell_count <- ncell(r)
  valid_count <- cell_count - na_count
  zcta_cell_count <- as.numeric(global(!is.na(zmask), "sum", na.rm = FALSE)[1, 1])
  zcta_na_count <- as.numeric(global(mask(is.na(r), zmask), "sum", na.rm = TRUE)[1, 1])
  zcta_valid_count <- zcta_cell_count - zcta_na_count
  vals <- global(r, c("min", "mean", "max"), na.rm = TRUE)

  png(
    filename = file.path(raw_map_dir, paste0(meta$pollutant, "_raw_grid.png")),
    width = 1600,
    height = 900,
    res = 160,
    bg = "white"
  )
  par(mar = c(3, 3, 4, 6))
  plot(
    r,
    col = viridis(100, option = "magma", direction = -1),
    colNA = "grey70",
    main = meta$label,
    axes = TRUE,
    plg = list(title = meta$units)
  )
  dev.off()

  summary_rows[[i]] <- data.table(
    pollutant = meta$pollutant,
    source_file = basename(meta$file),
    layer = meta$layer,
    cells = cell_count,
    valid_cells = valid_count,
    missing_cells = na_count,
    missing_fraction = na_count / cell_count,
    zcta_cells = zcta_cell_count,
    zcta_valid_cells = zcta_valid_count,
    zcta_missing_cells = zcta_na_count,
    zcta_missing_fraction = zcta_na_count / zcta_cell_count,
    min = vals[1, "min"],
    mean = vals[1, "mean"],
    max = vals[1, "max"]
  )
}

summary <- rbindlist(summary_rows)
fwrite(summary, file.path(raw_map_dir, "raw_grid_summary.csv"))
message("Raw grid maps written to ", normalizePath(raw_map_dir, winslash = "/"))
