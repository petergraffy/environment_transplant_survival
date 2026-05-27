library(ncdf4)
library(terra)

files <- c(
  pm = "C:/Users/Peter Graffy/Box/PM2.5/Monthly/Monthly/2005/V6GL02.04.CNNPM25.NA.200501-200501.nc",
  o3 = "C:/Users/Peter Graffy/Box/O3/CONUS_O3_MDA8_p01_monthly_mean_2005.nc",
  no2_old = "C:/Users/Peter Graffy/Box/NO2/Anenberg et al 2022 Global NO2/2005_final_1km.nc",
  no2_new = "C:/Users/Peter Graffy/Box/NO2/Anenberg et al 2022 Global NO2/annual_mean_tropomi_lur_conus_surface_no2_2024.v1.02.nc"
)

for (nm in names(files)) {
  cat("\n====", nm, "====\n")
  nc <- nc_open(files[[nm]])
  on.exit(nc_close(nc), add = TRUE)

  cat("dimensions:\n")
  for (d in nc$dim) {
    cat(
      " ", d$name,
      "len=", d$len,
      "units=", d$units,
      "first=", d$vals[1],
      "last=", d$vals[d$len],
      "\n"
    )
  }

  cat("variables:\n")
  for (v in nc$var) {
    cat(
      " ", v$name,
      "size=", paste(v$size, collapse = "x"),
      "units=", v$units,
      "dims=", paste(vapply(v$dim, function(x) x$name, character(1)), collapse = ","),
      "\n"
    )
  }
  nc_close(nc)

  r <- rast(files[[nm]])
  cat("terra names:", paste(names(r), collapse = ", "), "\n")
  cat("terra nlyr:", nlyr(r), "nrow:", nrow(r), "ncol:", ncol(r), "\n")
  cat("terra crs:", crs(r, describe = TRUE)$name, "\n")
  cat("terra ext:", as.character(ext(r)), "\n")
}
