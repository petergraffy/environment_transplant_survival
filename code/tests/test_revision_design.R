source("code/revision_2026/common.R")
stopifnot(identical(spell_ids(c(1, 2, 8, 12, 20), c(10, 3, 12, 15, 22)), c(1L,1L,1L,1L,2L)))
x <- concurrent_intervals(c(0,2,8), c(5,6,10), 0,10)
stopifnot(identical(x$concurrent_listings, c(1L,2L,1L,0L,1L)),
          sum(x$tstop-x$tstart) == 10, all(x$tstop > x$tstart))
stopifnot(identical(risk_counts(c(1,2,2,5), c(0,1,2,3,6)), c(4L,4L,3L,1L,0L)))
stopifnot(identical(spell_ids(c(1,5,11),c(5,10,12)),c(1L,1L,2L)))
stopifnot(nrow(concurrent_intervals(0,0,0,0))==0L)
clipped <- concurrent_intervals(c(-5,3,5),c(8,3,12),0,10)
stopifnot(sum(clipped$tstop-clipped$tstart)==10,
          identical(clipped$concurrent_listings,c(1L,2L,1L)))
future <- concurrent_intervals(c(0,5),c(10,10),0,10)
stopifnot(future$concurrent_listings[1]==1L,future$concurrent_listings[2]==2L)
cat("PASS: overlap chains, same-day boundaries, gaps, concurrency, and pre-event risk sets.\n")
