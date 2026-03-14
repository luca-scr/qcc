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

test_that("ProcessCapability returns expected structure", {
  chart <- qcc(make_grouped_xbar_data(), type = "xbar", nsigmas = 3)

  capability <- processCapability(chart, spec.limits = c(9.9, 10.1))

  expect_s3_class(capability, "processCapability")
  expect_equal(names(capability$spec.limits), c("LSL", "USL"))
  expect_equal(unname(capability$spec.limits), c(9.9, 10.1))

  expect_true(is.matrix(capability$indices))
  expect_equal(rownames(capability$indices), c("Cp", "Cp_l", "Cp_u", "Cp_k", "Cpm"))
  expect_equal(colnames(capability$indices)[1], "Value")
  expect_equal(ncol(capability$indices), 3)
})

test_that("ProcessCapability rejects invalid specs", {
  chart <- qcc(make_grouped_xbar_data(), type = "xbar", nsigmas = 3)

  expect_error(
    processCapability(chart, spec.limits = c(NA_real_, NaN)),
    "invalid specification limits"
  )
})

test_that("ProcessCapability returns expected indices", {
  chart <- qcc(make_grouped_xbar_data(), type = "xbar", nsigmas = 3)
  capability <- processCapability(chart, spec.limits = c(9.9, 10.1))

  expected_indices <- matrix(
    c(
      1.372667, 1.441300, 1.304033, 1.304033, 1.344463,
      0.808460, 0.911657, 0.820114, 0.727408, 0.804207,
      1.937713, 1.970943, 1.787953, 1.880659, 1.885321
    ),
    nrow = 5,
    ncol = 3,
    dimnames = list(
      c("Cp", "Cp_l", "Cp_u", "Cp_k", "Cpm"),
      c("Value", "2.5%", "97.5%")
    )
  )

  expect_equal(capability$indices, expected_indices, tolerance = 1e-6)
})

test_that("No visual regressions in process capability plots", {
  chart <- qcc(make_grouped_xbar_data(), type = "xbar", nsigmas = 3)
  capability <- processCapability(chart, spec.limits = c(9.9, 10.1))

  vdiffr::expect_doppelganger(
    "process capability plot",
    plot(capability)
  )
})

test_that("plot.processCapability fails with unexpected objects", {
  expect_error(
    plot.processCapability(make_grouped_xbar_data()),
    "an object of class `processCapability' is required"
  )
})
