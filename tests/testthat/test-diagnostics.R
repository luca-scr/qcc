test_that("qcc_diagnostics rejects non-qcc input", {
  expect_error(qcc_diagnostics(list()), "'object' must be of class 'qcc'")
  expect_error(qcc_diagnostics("not a qcc"), "'object' must be of class 'qcc'")
})

test_that("inactive patterns return NULL and do not compute", {
  q  <- qcc(rep(0, 20), type = "xbar.one", center = 0, std.dev = 1)
  dx <- qcc_diagnostics(q)
  expect_null(dx$pattern5)
  expect_null(dx$pattern6)
  expect_null(dx$pattern7)
  expect_null(dx$pattern8)
})

test_that("original qcc object is not modified", {
  q              <- qcc(c(-2, 2, rep(0, 15), 2, -2, 2),
                        type = "xbar.one", center = 0, std.dev = 1)
  orig_viol      <- q$violations
  orig_stats     <- q$statistics
  orig_center    <- q$center
  orig_std.dev   <- q$std.dev
  dx             <- qcc_diagnostics(q, pattern5 = TRUE, pattern8 = TRUE)
  expect_identical(q$violations, orig_viol)
  expect_identical(q$statistics, orig_stats)
  expect_identical(q$center,     orig_center)
  expect_identical(q$std.dev,    orig_std.dev)
  expect_identical(dx$qcc_object$violations, orig_viol)
})

test_that("result has class qcc_diagnostics", {
  q  <- qcc(rnorm(20), type = "xbar.one")
  dx <- qcc_diagnostics(q, pattern8 = TRUE)
  expect_s3_class(dx, "qcc_diagnostics")
})

# ---------------------------------------------------------------------------
# Pattern 5 — 15 consecutive in Zone C
# center=0, std.dev=1, nsigmas=3  →  Zone C = [-1, 1]
# Series: 2 outside points, 15 in Zone C, 3 outside points
# First all-in-Zone-C window ends at position 17
# ---------------------------------------------------------------------------

test_that("pattern5 detects 15 consecutive in Zone C", {
  s  <- c(-2, 2, rep(0, 15), 2, -2, 2)
  q  <- qcc(s, type = "xbar.one", center = 0, std.dev = 1)
  dx <- qcc_diagnostics(q, pattern5 = TRUE)
  expect_type(dx$pattern5, "integer")
  expect_equal(dx$pattern5, 17L)
})

test_that("pattern5 returns integer(0) when no run of 15 in Zone C", {
  # All points outside Zone C
  s  <- rep(c(2, -2), 10)
  q  <- qcc(s, type = "xbar.one", center = 0, std.dev = 1)
  dx <- qcc_diagnostics(q, pattern5 = TRUE)
  expect_equal(dx$pattern5, integer(0L))
})

test_that("pattern5 returns integer(0) for series shorter than 15", {
  s  <- rep(0, 14)
  q  <- qcc(s, type = "xbar.one", center = 0, std.dev = 1)
  dx <- qcc_diagnostics(q, pattern5 = TRUE)
  expect_equal(dx$pattern5, integer(0L))
})

# ---------------------------------------------------------------------------
# Pattern 6 — 8 consecutive outside Zone C
# Series: 5 in Zone C, 8 outside, 5 in Zone C
# First all-outside window ends at position 13
# ---------------------------------------------------------------------------

test_that("pattern6 detects 8 consecutive outside Zone C", {
  s  <- c(rep(0, 5), rep(2, 8), rep(0, 5))
  q  <- qcc(s, type = "xbar.one", center = 0, std.dev = 1)
  dx <- qcc_diagnostics(q, pattern6 = TRUE)
  expect_type(dx$pattern6, "integer")
  expect_equal(dx$pattern6, 13L)
})

test_that("pattern6 detects run on negative side", {
  s  <- c(rep(0, 5), rep(-2, 8), rep(0, 5))
  q  <- qcc(s, type = "xbar.one", center = 0, std.dev = 1)
  dx <- qcc_diagnostics(q, pattern6 = TRUE)
  expect_equal(dx$pattern6, 13L)
})

test_that("pattern6 returns integer(0) when no run of 8 outside Zone C", {
  s  <- rep(0, 20)
  q  <- qcc(s, type = "xbar.one", center = 0, std.dev = 1)
  dx <- qcc_diagnostics(q, pattern6 = TRUE)
  expect_equal(dx$pattern6, integer(0L))
})

test_that("pattern6 returns integer(0) for series shorter than 8", {
  s  <- rep(2, 7)
  q  <- qcc(s, type = "xbar.one", center = 0, std.dev = 1)
  dx <- qcc_diagnostics(q, pattern6 = TRUE)
  expect_equal(dx$pattern6, integer(0L))
})

# ---------------------------------------------------------------------------
# Pattern 7 — 14 consecutive alternating points
# Series: 5 zeros + 14 alternating (3,-3 repeated) + 5 zeros
# Alternating region at positions 6-19; detections within that range
# ---------------------------------------------------------------------------

test_that("pattern7 detects 14 consecutive alternating points", {
  s  <- c(rep(0, 5), rep(c(3, -3), 7), rep(0, 5))
  q  <- qcc(s, type = "xbar.one", center = 0, std.dev = 1)
  dx <- qcc_diagnostics(q, pattern7 = TRUE)
  expect_type(dx$pattern7, "integer")
  expect_gt(length(dx$pattern7), 0L)
  expect_true(all(dx$pattern7 >= 18L & dx$pattern7 <= 20L))
})

test_that("pattern7 returns integer(0) for monotone series", {
  # Strictly increasing: no alternation at all
  s  <- 1:25
  q  <- qcc(s, type = "xbar.one", center = 13, std.dev = 5)
  dx <- qcc_diagnostics(q, pattern7 = TRUE)
  expect_equal(dx$pattern7, integer(0L))
})

test_that("pattern7 returns integer(0) for series shorter than 14", {
  s  <- rep(c(1, -1), 6)   # 12 points
  q  <- qcc(s, type = "xbar.one", center = 0, std.dev = 1)
  dx <- qcc_diagnostics(q, pattern7 = TRUE)
  expect_equal(dx$pattern7, integer(0L))
})

# ---------------------------------------------------------------------------
# Pattern 8 — 7 consecutive trending points
# Series: 5 zeros + 7 strictly increasing + 5 zeros  (17 points)
# diffs[5:10] all positive → window ends at index 11
# ---------------------------------------------------------------------------

test_that("pattern8 detects 7 consecutive upward trend", {
  s  <- c(rep(0, 5), 1:7, rep(0, 5))
  q  <- qcc(s, type = "xbar.one", center = 0, std.dev = 1)
  dx <- qcc_diagnostics(q, pattern8 = TRUE)
  expect_type(dx$pattern8, "integer")
  expect_true(11L %in% dx$pattern8)
})

test_that("pattern8 detects 7 consecutive downward trend", {
  s  <- c(rep(0, 5), 7:1, rep(0, 5))
  q  <- qcc(s, type = "xbar.one", center = 0, std.dev = 1)
  dx <- qcc_diagnostics(q, pattern8 = TRUE)
  expect_true(11L %in% dx$pattern8)
})

test_that("pattern8 returns integer(0) for flat series", {
  s  <- rep(0, 20)
  q  <- qcc(s, type = "xbar.one", center = 0, std.dev = 1)
  dx <- qcc_diagnostics(q, pattern8 = TRUE)
  expect_equal(dx$pattern8, integer(0L))
})

test_that("pattern8 returns integer(0) for series shorter than 7", {
  s  <- 1:6
  q  <- qcc(s, type = "xbar.one", center = 3.5, std.dev = 2)
  dx <- qcc_diagnostics(q, pattern8 = TRUE)
  expect_equal(dx$pattern8, integer(0L))
})

test_that("pattern8 does not trigger on trend broken by a tie", {
  # 6 increasing then a repeat (tie breaks the run)
  s  <- c(rep(0, 5), 1, 2, 3, 4, 5, 5, 6, rep(0, 5))
  q  <- qcc(s, type = "xbar.one", center = 0, std.dev = 1)
  dx <- qcc_diagnostics(q, pattern8 = TRUE)
  # The run 1,2,3,4,5 is only 5 points; the run after the tie has fewer than 7
  expect_equal(dx$pattern8, integer(0L))
})

# ---------------------------------------------------------------------------
# Zone C infrastructure
# ---------------------------------------------------------------------------

test_that("Zone C limits are strictly inside control limits", {
  q      <- qcc(rnorm(30), type = "xbar.one")
  dx     <- qcc_diagnostics(q, pattern5 = TRUE)
  zone_c <- dx$zone_c
  lcl    <- q$limits[, 1]
  ucl    <- q$limits[, 2]
  expect_true(all(zone_c[, 1] > lcl))
  expect_true(all(zone_c[, 2] < ucl))
})

test_that("Zone C is symmetric around center", {
  q      <- qcc(rnorm(30), type = "xbar.one")
  dx     <- qcc_diagnostics(q, pattern5 = TRUE)
  zone_c <- dx$zone_c
  center <- q$center
  expect_equal(zone_c[, 2] - center, center - zone_c[, 1], tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Phase II (newdata) support
# ---------------------------------------------------------------------------

test_that("stats field concatenates phase I and phase II", {
  s1 <- rnorm(20)
  s2 <- rnorm(10)
  q  <- qcc(s1, type = "xbar.one", newdata = s2)
  dx <- qcc_diagnostics(q, pattern8 = TRUE)
  expect_equal(length(dx$stats), 30L)
})

# ---------------------------------------------------------------------------
# Multiple patterns simultaneously
# ---------------------------------------------------------------------------

test_that("multiple patterns can be requested in a single call", {
  s  <- c(-2, 2, rep(0, 15), 2, -2, 2)
  q  <- qcc(s, type = "xbar.one", center = 0, std.dev = 1)
  dx <- qcc_diagnostics(q, pattern5 = TRUE, pattern6 = TRUE,
                          pattern7 = TRUE, pattern8 = TRUE)
  expect_type(dx$pattern5, "integer")
  expect_type(dx$pattern6, "integer")
  expect_type(dx$pattern7, "integer")
  expect_type(dx$pattern8, "integer")
})
