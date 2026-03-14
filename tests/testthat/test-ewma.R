test_that("ewma computes smoothing, sigma, and limits from explicit parameters", {
  data <- matrix(c(11, 13, 12), ncol = 1)
  expected_sigma <- c(1, sqrt(1.25), sqrt(1.3125))
  expected_limits <- cbind(
    LCL = 10 - 2 * expected_sigma,
    UCL = 10 + 2 * expected_sigma
  )

  chart <- ewma(
    data,
    sizes = 1,
    center = 10,
    std.dev = 2,
    lambda = 0.5,
    nsigmas = 2
  )

  expect_equal(unname(chart$ewma), c(10.5, 11.75, 11.875), tolerance = 1e-12)
  expect_equal(unname(chart$sigma), expected_sigma, tolerance = 1e-12)
  expect_equal(unname(chart$limits), unname(expected_limits), tolerance = 1e-12)
})

test_that("ewmaSmooth orders x and smooths y in sorted order", {
  smoothed <- ewmaSmooth(
    x = c(3, 1, 2),
    y = c(30, 10, 20),
    lambda = 0.5,
    start = 0
  )

  expect_equal(smoothed$x, c(1, 2, 3))
  expect_equal(smoothed$y, c(5, 12.5, 21.25), tolerance = 1e-12)
})

test_that("ewma validates lambda range", {
  data <- matrix(c(1, 2, 3), ncol = 1)

  expect_error(
    ewma(data, sizes = 1, center = 0, std.dev = 1, lambda = 1.2),
    "lambda parameter must be between 0 and 1"
  )

  expect_error(
    ewma(data, sizes = 1, center = 0, std.dev = 1, lambda = -0.2),
    "lambda parameter must be between 0 and 1"
  )
})

test_that("ewma phase-II extension propagates labels and flags violations", {
  phase_i <- matrix(c(0, 0), ncol = 1)
  phase_ii <- matrix(c(0, 5), ncol = 1, dimnames = list(c("trial-A", "trial-B"), NULL))

  chart <- ewma(
    phase_i,
    sizes = 1,
    center = 0,
    std.dev = 1,
    lambda = 1,
    nsigmas = 2,
    newdata = phase_ii,
    newsizes = 1
  )

  expect_equal(names(chart$newstats), c("trial-A", "trial-B"))
  expect_equal(names(chart$ewma), c("1", "2", "trial-A", "trial-B"))
  expect_equal(unname(chart$ewma), c(0, 0, 0, 5), tolerance = 1e-12)
  expect_equal(unname(chart$violations), c(NA_real_, NA_real_, NA_real_, 1))
})

test_that("no visual regressions ewma chart plots", {
  data <- matrix(c(11, 13, 12), ncol = 1)
  chart <- ewma(
    data,
    sizes = 1,
    center = 10,
    std.dev = 2,
    lambda = 0.5,
    nsigmas = 2
  )

  vdiffr::expect_doppelganger(
    "ewma plot - no fill",
    plot(chart, fill = FALSE)
  )
  vdiffr::expect_doppelganger(
    "ewma plot - chartall",
    plot(chart, chart.all = TRUE)
  )
})

test_that("no visual regressions in phase II ewma chart plots", {
  phase_i <- matrix(c(0, 0), ncol = 1)
  phase_ii <- matrix(
    c(0, 5),
    ncol = 1,
    dimnames = list(c("trial-A", "trial-B"), NULL)
  )

  chart <- ewma(
    phase_i,
    sizes = 1,
    center = 0,
    std.dev = 1,
    lambda = 1,
    nsigmas = 2,
    newdata = phase_ii,
    newsizes = 1
  )

  vdiffr::expect_doppelganger(
    "ewma phase II plot",
    plot(chart)
  )
})
