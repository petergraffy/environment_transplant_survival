source("code/daily_pollution_inputs.R")

scripts <- c(list.files("code", pattern = "^(78|79|80|81|84|86|88|90).*R$", full.names = TRUE),
             "code/figures_tables/60_plot_primary_pollution_study_period_maps.R",
             "code/figures_tables/92_plot_pollution_absolute_1yr_survival_maps.R",
             "code/prepare_daily_pollution_inputs.R")
invisible(lapply(scripts, parse))

local({
  root <- tempfile()
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  file <- file.path(root, "source.parquet")
  writeBin(as.raw(1:3), file)
  key <- function() daily_pollution_cache_path("exposure.csv.gz", file, "value", "output", "annual")
  first <- key()
  stopifnot(identical(first, key()))
  writeBin(as.raw(1:4), file)
  stopifnot(!identical(first, key()))
  missing_path <- file.path(root, "o3_zcta_daily_parquet")
  dir.create(missing_path)
  stopifnot(inherits(try(daily_pollution_files(missing_path), silent = TRUE), "try-error"))
})

for (directory in daily_pollution_releases$directory) {
  files <- daily_pollution_files(file.path("data", "release", directory))
  stopifnot(length(files) == 20L)
}
stopifnot(daily_pollution_end_date == as.Date("2024-12-31"))
cat("PASS: scripts parse; cache keys are stable and change with inputs; incomplete releases fail; coverage ends in 2024.\n")
