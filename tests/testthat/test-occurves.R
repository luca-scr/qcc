make_grouped_xbar_data <- function() {
  matrix(
    c(
      10.00, 10.02, 9.98, 10.01,
      10.03, 9.99, 10.01, 10.00,
      9.97, 10.04, 10.02, 9.99
    ),
    nrow = 3,
    byrow = TRUE
  )
}

test_that("ocCurves on xbar qcc object returns coherent beta and ARL matrices", {
  chart <- qcc(make_grouped_xbar_data(), type = "xbar", nsigmas = 3)

  oc <- ocCurves(chart)

  expect_s3_class(oc, "ocCurves")
  expect_equal(oc$type, "xbar")
  expect_true(is.matrix(oc$beta))
  expect_true(is.matrix(oc$ARL))

  expect_equal(nrow(oc$beta), length(oc$shift))
  expect_equal(ncol(oc$beta), length(oc$size))
  expect_equal(dim(oc$ARL), dim(oc$beta))
  expect_equal(oc$ARL, 1 / (1 - oc$beta), tolerance = 1e-12)
})

test_that("ocCurves rejects unsupported qcc chart type", {
  one_at_time <- matrix(c(10.1, 9.9, 10.0, 10.2, 9.8), ncol = 1)
  chart <- qcc(one_at_time, type = "xbar.one")

  expect_error(
    ocCurves(chart),
    "Operating characteristic curves not available for this type of chart"
  )
})
