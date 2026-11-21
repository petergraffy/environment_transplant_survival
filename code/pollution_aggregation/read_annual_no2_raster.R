# Use the NetCDF coordinates to orient annual TROPOMI values. GDAL can
# reverse the latitude rows in v2 files while keeping the correct extent.
read_annual_no2_raster <- function(path) {
  r <- terra::rast(path)
  if (!identical(names(r), "surface_no2")) return(r)

  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc))
  variable <- nc$var$surface_no2
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  lat_order <- order(lat, decreasing = TRUE)
  stopifnot(identical(vapply(variable$dim, function(d) d$name, ""), c("lon", "lat")),
    length(lon) == ncol(r), length(lat) == nrow(r),
    all(diff(lon) > 0), !anyDuplicated(lat),
    max(abs(terra::xFromCol(r, seq_len(ncol(r))) - lon)) < 1e-5,
    max(abs(terra::yFromRow(r, seq_len(nrow(r))) - lat[lat_order])) < 1e-5)
  terra::values(r) <- as.vector(ncdf4::ncvar_get(nc, "surface_no2")[, lat_order, drop = FALSE])
  r
}
