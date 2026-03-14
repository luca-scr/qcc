test_that("cusum recursion and violation checks use strict decision boundaries", {
  data <- matrix(c(1.5, 0.5, 0.6, -1.5, -0.5, -0.8), ncol = 1)

  chart <- cusum(
    data,
    sizes = 1,
    center = 0,
    std.dev = 1,
    decision.interval = 1,
    se.shift = 1,
    head.start = 0
  )

  expect_equal(unname(chart$pos), c(1.0, 1.0, 1.1, 0.0, 0.0, 0.0), tolerance = 1e-12)
  expect_equal(unname(chart$neg), c(0.0, 0.0, 0.0, -1.0, -1.0, -1.3), tolerance = 1e-12)
  expect_equal(
    chart$violations$upper,
    c(NA_real_, NA_real_, 1, NA_real_, NA_real_, NA_real_)
  )
  expect_equal(
    chart$violations$lower,
    c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, 1)
  )
})

test_that("cusum phase-II data are appended with sequential labels", {
  phase_i <- matrix(c(10, 12), ncol = 1)
  phase_ii <- matrix(c(14, 16), ncol = 1)

  chart <- cusum(
    phase_i,
    sizes = 1,
    center = 10,
    std.dev = 2,
    decision.interval = 4,
    se.shift = 1,
    newdata = phase_ii,
    newsizes = 1
  )

  expect_equal(unname(chart$newstats), c(14, 16))
  expect_equal(names(chart$newstats), c("3", "4"))
  expect_equal(length(chart$pos), 4)
})

test_that("cusum phase-II keeps rowname labels and uses explicit newsizes", {
  phase_i <- matrix(c(10, 10, 10, 10, 10, 10), nrow = 2, byrow = TRUE)
  phase_ii <- matrix(
    c(14, 14, NA, 16, 16, 16),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("lot-A", "lot-B"), NULL)
  )

  chart <- cusum(
    phase_i,
    sizes = 3,
    center = 10,
    std.dev = 2,
    decision.interval = 10,
    se.shift = 1,
    newdata = phase_ii,
    newsizes = c(2, 3)
  )

  expect_equal(names(chart$newstats), c("lot-A", "lot-B"))
  expect_equal(chart$newsizes, c(2, 3))
  expect_equal(unname(chart$pos), c(0, 0, 2.32842712474619, 7.02457954745282), tolerance = 1e-12)
})

test_that("cusum validates head.start bounds", {
  data <- matrix(c(1, 2), ncol = 1)

  expect_error(
    cusum(
      data,
      sizes = 1,
      center = 0,
      std.dev = 1,
      decision.interval = 2,
      head.start = 2
    ),
    "head.start must be non-negative and less than decision.interval"
  )
})

test_that("no visual regressions in cusum plots", {
    phase_i <- matrix(c(10, 12), ncol = 1)
  phase_ii <- matrix(c(14, 16), ncol = 1)

  chart <- cusum(
    phase_i,
    sizes = 1,
    center = 10,
    std.dev = 2,
    decision.interval = 4,
    se.shift = 1,
    newdata = phase_ii,
    newsizes = 1
  )

  vdiffr::expect_doppelganger(
    "cusum phase II plot",
    plot(chart)
  )
  succeed()

})
