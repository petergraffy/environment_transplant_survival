source("code/r_runtime.R"); ensure_user_library()
suppressPackageStartupMessages({library(data.table); library(survival)})
source("code/lung_survival_benefit_helpers.R")
# Closed-form area for a step hazard, including a shifted decision date and an
# event exactly at the restriction horizon (which contributes no further area).
h <- data.table(strata = c("a", "a", "b"), time = c(10, 365, 20), hazard = c(.2, .4, .3))
rr <- matrix(c(1, 2, .5, 3), nrow = 2)
got <- lung_rmst_integrate(h, c("a", "b"), c(5, 0), rr, c(0, 10, 365))
expected <- sapply(c(0, 10, 365), function(t) {
  a <- pmin(t, 5) + pmax(pmin(t, 360) - 5, 0) * exp(-.2 * rr[1, ]) + pmax(t - 360, 0) * exp(-.4 * rr[1, ])
  b <- pmin(t, 20) + pmax(t - 20, 0) * exp(-.3 * rr[2, ])
  (a + b) / 2
})
stopifnot(max(abs(got - expected)) < 1e-10)
shifted <- lung_rmst_integrate(h, "a", 15, matrix(2, nrow = 1), 365)
stopifnot(abs(shifted[1, 1] - (350 + 15 * exp(-.2 * 2))) < 1e-10)
# No hazard jumps implies exactly the available horizon for every scenario.
none <- data.table(strata = "a", time = 1000, hazard = .1)
stopifnot(all(lung_rmst_integrate(none, "a", 0, matrix(c(1, 5), nrow = 1), c(0, 365)) ==
  matrix(c(0, 0, 365, 365), nrow = 2)))
# Identical fitted pre/post models, age definitions and origins imply zero RMST
# benefit; this also checks concentration/horizon ordering and covariate prediction.
set.seed(93)
n <- 200L
d <- data.table(PERS_ID = seq_len(n), tstart = 0, tstop = runif(n, 10, 1000), event = sample(0:1, n, TRUE),
  pm25 = runif(n, 4, 12), age = runif(n, 20, 70), sex = sample(c("F", "M"), n, TRUE),
  race = sample(c("A", "B"), n, TRUE), zcta_svi_proxy = runif(n), listing_center = sample(c("a", "b"), n, TRUE),
  las_per10 = runif(n, 3, 6), las_unavailable = 0L)
m <- lung_fit_cox(d, "pm25")
r <- copy(d[1:5]); r[, `:=`(waiting_days = 0, listing_age = age, transplant_age = age)]
q <- lung_predict_rmst(m, m, r, "pm25", c(5, 10), c(0, 180, 365))
stopifnot(nrow(q) == 6L, all(abs(q$rmst_benefit_days) < 1e-10),
  identical(q$concentration, rep(c(5, 10), 3)), identical(q$horizon_days, rep(c(0, 180, 365), each = 2)))
cat("RMST tests passed: exact step integration, shifted origins, horizon knots, zero hazard and paired-model predictions.\n")
