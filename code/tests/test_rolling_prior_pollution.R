source("code/r_runtime.R")
ensure_user_library()
suppressPackageStartupMessages(library(data.table))
source("code/rolling_prior_pollution.R")

dates <- seq(as.Date("2019-01-01"), as.Date("2021-01-02"), by = "day")
d <- data.frame(zip = "01234", date = dates, value = seq_along(dates))
q <- data.frame(zip = c("01234", "01234", "99999"),
                 start = as.Date(c("2019-07-02", "2020-01-02", "2020-01-02")),
                 end = as.Date(c("2020-06-30", "2020-12-31", "2020-12-31")))
s <- rolling_window_sum(d, q)
for (i in 1:2) {
  brute <- d$value[d$date >= q$start[i] & d$date <= q$end[i]]
  stopifnot(s$count[i] == 365, s$total[i] == sum(brute))
}
stopifnot(s$count[3] == 0, s$total[3] == 0)
# Partitioning at New Year must give the same result as reading both years.
s1 <- rolling_window_sum(d[d$date < as.Date("2020-01-01"), ], q)
s2 <- rolling_window_sum(d[d$date >= as.Date("2020-01-01"), ], q)
stopifnot(all(s1$total + s2$total == s$total), all(s1$count + s2$count == s$count))
d$value[d$date == as.Date("2020-02-29")] <- NA_real_
stopifnot(rolling_window_sum(d, q)$count[2] == 364)
stopifnot(inherits(try(rolling_window_sum(rbind(d, d[1, ]), q), silent = TRUE), "try-error"))

m <- data.frame(zip = "01234", year = 2019L, month = 1:12, no2_ppbv = 1:12)
a <- data.frame(zip = "01234", year = 2018L, no2_annual = 10)
nq <- data.frame(zip = "01234", index_date = as.Date(c("2020-01-01", "2020-01-31", "2019-07-15")))
n <- rolling_no2_windows(nq, m, a)
month_days <- c(31,28,31,30,31,30,31,31,30,31,30,31)
stopifnot(all.equal(n$no2_prior_ppb[1], weighted.mean(1:12, month_days)),
           n$no2_prior_ppb[1] == n$no2_prior_ppb[2],
           all(n$no2_prior_end < nq$index_date),
           n$no2_prior_annual_months[3] == 6,
           n$no2_prior_start[3] == as.Date("2018-07-01"),
           n$no2_prior_end[3] == as.Date("2019-06-30"))
stopifnot(is.na(rolling_no2_windows(nq, m, annual_fallback = FALSE)$no2_prior_ppb[3]))
# A hole within observed monthly coverage must not be filled using an annual value.
stopifnot(is.na(rolling_no2_windows(nq[1, ], m[-5, ], rbind(a, data.frame(zip="01234",year=2019L,no2_annual=99)))$no2_prior_ppb))
paths <- c(list.files("code", pattern = "^(78|80|81|86|88|90).*R$", full.names = TRUE),
           "code/figures_tables/92_plot_pollution_absolute_1yr_survival_maps.R")
invisible(lapply(paths, parse))
cat("PASS: rolling daily windows, leap days, missing days, duplicate dates, year partitions, 12-month windows, annual fallback, and script parsing.\n")
