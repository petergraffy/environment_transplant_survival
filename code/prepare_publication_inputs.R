#!/usr/bin/env Rscript
source("code/r_runtime.R")
ensure_user_library()

inputs <- read.csv("data/reference/publication_inputs.csv", colClasses = "character")
download <- "--download" %in% commandArgs(trailingOnly = TRUE)
options(timeout = max(600, getOption("timeout")))
verify <- function(path, size, hash) {
  file.exists(path) && file.size(path) == as.numeric(size) && unname(tools::sha256sum(path)) == hash
}
for (i in seq_len(nrow(inputs))) {
  x <- inputs[i, ]
  if (verify(x$path, x$size_bytes, x$sha256)) next
  if (file.exists(x$path)) stop("Input differs from the recorded checksum: ", x$path)
  if (!download) stop("Missing input: ", x$path, ". Run the inputs stage with --download.")
  dir.create(dirname(x$path), recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(tmpdir = dirname(x$path))
  message("Downloading ", basename(x$path))
  tryCatch({
    status <- download.file(x$url, stage, mode = "wb", quiet = TRUE)
    if (status != 0L || !verify(stage, x$size_bytes, x$sha256)) stop("Invalid download for ", x$path)
    if (!file.rename(stage, x$path)) stop("Could not finalize ", x$path)
  }, finally = unlink(stage))
}

boundary_dir <- "data/cache/cb_2020_us_zcta520_500k"
boundary_files <- file.path(boundary_dir, paste0("cb_2020_us_zcta520_500k.", c("shp", "shx", "dbf", "prj")))
if (!all(file.exists(boundary_files))) {
  unzip(paste0(boundary_dir, ".zip"), exdir = boundary_dir)
  stopifnot(all(file.exists(boundary_files)))
}
message("Verified ", nrow(inputs), " public input files.")
