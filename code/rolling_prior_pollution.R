# Listing-date-specific baseline exposures. Time-varying inputs are unchanged.
source(file.path("code", "daily_pollution_inputs.R"))

rolling_window_sum <- function(daily, queries) {
  # Inclusive endpoints, exact observed-day counts, and no listing-day exposure.
  daily <- data.table::as.data.table(daily)
  daily <- data.table::copy(daily)
  daily[, date := as.Date(date)]
  data.table::setorder(daily, zip, date)
  if (anyDuplicated(daily[, .(zip, date)])) stop("Duplicate ZCTA-date exposure records")
  daily[, `:=`(total = cumsum(data.table::fifelse(is.finite(value), value, 0)),
                count = cumsum(as.integer(is.finite(value)))), by = zip]
  data.table::setkey(daily, zip, date)
  q <- data.table::as.data.table(queries)
  end <- daily[q[, .(zip, date = end)], on = .(zip, date), roll = Inf]
  before <- daily[q[, .(zip, date = start - 1)], on = .(zip, date), roll = Inf]
  data.frame(total = data.table::fcoalesce(end$total, 0) - data.table::fcoalesce(before$total, 0),
             count = data.table::fcoalesce(end$count, 0L) - data.table::fcoalesce(before$count, 0L))
}

rolling_input_key <- function(queries, files, mode) {
  tmp <- tempfile()
  on.exit(unlink(tmp))
  info <- file.info(files)
  saveRDS(list(version = 1L, queries = as.data.frame(queries),
               files = normalizePath(files, winslash = "/"), size = info$size,
               mtime = as.numeric(info$mtime), mode = mode), tmp)
  unname(tools::md5sum(tmp))
}

read_rolling_daily <- function(queries, pollutant) {
  spec <- daily_pollution_releases[daily_pollution_releases$pollutant == pollutant, ]
  files <- daily_pollution_files(file.path("data", "release", spec$directory))
  cache <- file.path("data", "cache", "rolling_prior_pollution",
                     paste0(pollutant, "_", rolling_input_key(queries, files, "365_complete_days"), ".rds"))
  if (file.exists(cache)) return(readRDS(cache))
  q <- data.table::as.data.table(data.table::copy(queries))
  q[, `:=`(start = index_date - 365, end = index_date - 1, total = 0, count = 0L)]
  for (file in files) {
    year <- as.integer(sub(".*_([0-9]{4})[.]parquet$", "\\1", file))
    first <- as.Date(paste0(year, "-01-01"))
    last <- as.Date(paste0(year, "-12-31"))
    ids <- which(!is.na(q$index_date) & !is.na(q$zip) & q$start <= last & q$end >= first)
    if (!length(ids)) next
    message("Rolling ", pollutant, ": reading ", year, " for ", length(ids), " ZCTA/listing-date windows")
    daily <- data.table::as.data.table(arrow::read_parquet(file, col_select = c("zip", "date", spec$value_column)))
    data.table::setnames(daily, spec$value_column, "value")
    daily <- daily[zip %in% q$zip[ids]]
    sums <- rolling_window_sum(daily, q[ids, .(zip, start = pmax(start, first), end = pmin(end, last))])
    q[ids, `:=`(total = total + sums$total, count = count + sums$count)]
    rm(daily, sums)
    gc(FALSE)
  }
  q[, value := data.table::fifelse(count == 365L, total / 365, NA_real_)]
  result <- q[, .(zip, index_date, value, covered_days = count,
                   window_start = start, window_end = end)]
  data.table::setnames(result, c("value", "covered_days", "window_start", "window_end"),
                       c(paste0(pollutant, "_prior_", if (pollutant == "pm25") "ug_m3" else "ppb"),
                         paste0(pollutant, c("_prior_days", "_prior_start", "_prior_end"))))
  dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
  saveRDS(as.data.frame(result), cache)
  as.data.frame(result)
}

rolling_no2_windows <- function(queries, monthly, annual = NULL, annual_fallback = TRUE) {
  q <- data.table::as.data.table(data.table::copy(queries))
  q[, id := .I]
  q[, listing_month := as.Date(format(index_date, "%Y-%m-01"))]
  q[, month_id := as.integer(format(listing_month, "%Y")) * 12L + as.integer(format(listing_month, "%m")) - 1L]
  windows <- q[rep(seq_len(.N), each = 12L), .(id, zip, index_date, month_id)]
  windows[, month_id := month_id - rep(12:1, times = nrow(q))]
  windows[, `:=`(year = month_id %/% 12L, month = month_id %% 12L + 1L)]
  windows[, start := as.Date(sprintf("%04d-%02d-01", year, month))]
  windows[, end := as.Date(sprintf("%04d-%02d-01", (month_id + 1L) %/% 12L,
                                   (month_id + 1L) %% 12L + 1L)) - 1]
  windows[, days := as.integer(end - start) + 1L]
  monthly <- data.table::as.data.table(monthly)
  if (anyDuplicated(monthly[, .(zip, year, month)])) stop("Duplicate NO2 ZCTA-month records")
  first_month <- min(monthly$year * 12L + monthly$month - 1L)
  windows <- monthly[windows, on = .(zip, year, month)]
  windows[, from_annual := FALSE]
  if (annual_fallback) {
    annual <- data.table::as.data.table(annual)
    if (anyDuplicated(annual[, .(zip, year)])) stop("Duplicate annual NO2 records")
    windows <- annual[windows, on = .(zip, year)]
    # Annual values approximate months only before the monthly source begins.
    windows[, from_annual := month_id < first_month & is.finite(no2_annual)]
    windows[from_annual == TRUE, no2_ppbv := no2_annual]
  }
  windows[, .(
    zip = zip[1], index_date = index_date[1],
    no2_prior_ppb = if (.N == 12L && all(is.finite(no2_ppbv))) sum(no2_ppbv * days) / sum(days) else NA_real_,
    no2_prior_months = sum(is.finite(no2_ppbv)),
    no2_prior_annual_months = sum(from_annual),
    no2_prior_start = min(start), no2_prior_end = max(end)
  ), by = id][, id := NULL][] |> as.data.frame()
}

read_rolling_prior_pollution <- function(dat, annual_fallback = TRUE) {
  queries <- unique(data.frame(zip = sprintf("%05d", as.integer(dat$candidate_zip)),
                               index_date = as.Date(dat$index_date)))
  queries <- queries[!is.na(suppressWarnings(as.integer(queries$zip))) & !is.na(queries$index_date), ]
  queries <- queries[order(queries$zip, queries$index_date), ]
  rownames(queries) <- NULL
  if (!nrow(queries)) stop("No valid ZCTA/listing-date pairs")
  monthly_files <- sort(list.files("data/release/no2_zcta_monthly_parquet", "[.]parquet$", full.names = TRUE))
  annual_file <- "data/release/air_pollution_zcta_parquet/air_pollution_zcta_no2_annual_2005_2025.parquet"
  if (!length(monthly_files)) stop("Monthly NO2 release is missing")
  key <- rolling_input_key(queries, c(monthly_files, if (annual_fallback) annual_file),
                           paste0("12_complete_day_weighted_months_annual_", annual_fallback))
  cache <- file.path("data", "cache", "rolling_prior_pollution", paste0("no2_", key, ".rds"))
  if (file.exists(cache)) {
    no2 <- readRDS(cache)
  } else {
    monthly <- dplyr::collect(dplyr::select(arrow::open_dataset(monthly_files), zip, year, month, no2_ppbv))
    annual <- if (annual_fallback) {
      dplyr::transmute(arrow::read_parquet(annual_file), zip, year, no2_annual = no2)
    } else NULL
    chunks <- split(seq_len(nrow(queries)), ceiling(seq_len(nrow(queries)) / 25000L))
    no2 <- dplyr::bind_rows(lapply(chunks, function(ids) {
      message("Rolling NO2: windows ", min(ids), "-", max(ids))
      rolling_no2_windows(queries[ids, ], monthly, annual, annual_fallback)
    }))
    dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
    saveRDS(no2, cache)
  }
  list(pm25 = read_rolling_daily(queries, "pm25"), o3 = read_rolling_daily(queries, "o3"), no2 = no2)
}
