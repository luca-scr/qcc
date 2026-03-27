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

# TODO: ProcessCapability returns expected indices. One-sided case
# TODO: processCapability returns expected indices. target-off mean
test_that("ProcessCapability returns expected indices. Two-sided case", {
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

test_that("ProcessCapability returns valid Cp-family indices for one- and two-sided specs, and Cpm respects the target", {
  chart <- qcc(make_grouped_xbar_data(), type = "xbar", nsigmas = 3)

  two_sided_specs <- c(9.9, 10.1)
  off_target <- 10.02
  midpoint_target <- mean(two_sided_specs)

  get_indices <- function(spec.limits, target) {
    processCapability(
      chart,
      spec.limits = spec.limits,
      target = target
    )$indices[, "Value"]
  }

  upper_indices <- get_indices(c(NA_real_, two_sided_specs[2]), off_target)

  expect_true(is.na(upper_indices["Cp"]))
  expect_true(is.na(upper_indices["Cp_l"]))
  expect_false(is.na(upper_indices["Cp_u"]))
  expect_false(is.na(upper_indices["Cp_k"]))
  expect_equal(upper_indices[["Cp_k"]], upper_indices[["Cp_u"]], tolerance = 1e-6)
  expect_true(is.na(upper_indices["Cpm"]))

  lower_indices <- get_indices(c(two_sided_specs[1], NA_real_), off_target)

  expect_true(is.na(lower_indices["Cp"]))
  expect_true(is.na(lower_indices["Cp_u"]))
  expect_false(is.na(lower_indices["Cp_l"]))
  expect_false(is.na(lower_indices["Cp_k"]))
  expect_equal(lower_indices[["Cp_k"]], lower_indices[["Cp_l"]], tolerance = 1e-6)
  expect_true(is.na(lower_indices["Cpm"]))

  two_sided_off_indices <- get_indices(two_sided_specs, off_target)
  two_sided_mid_indices <- get_indices(two_sided_specs, midpoint_target)

  cp_family <- c("Cp", "Cp_l", "Cp_u", "Cp_k")

  expect_false(anyNA(two_sided_off_indices[cp_family]))
  expect_false(anyNA(two_sided_mid_indices[cp_family]))

  expect_equal(
    two_sided_off_indices[["Cp_k"]],
    min(two_sided_off_indices[["Cp_l"]], two_sided_off_indices[["Cp_u"]]),
    tolerance = 1e-6
  )

  expect_equal(
    two_sided_mid_indices[["Cp_k"]],
    min(two_sided_mid_indices[["Cp_l"]], two_sided_mid_indices[["Cp_u"]]),
    tolerance = 1e-6
  )

  # Cp-family indices should not depend on target
  expect_equal(
    unname(two_sided_off_indices[cp_family]),
    unname(two_sided_mid_indices[cp_family]),
    tolerance = 1e-6
  )

  # Cpm should remain finite
  expect_false(is.na(two_sided_off_indices["Cpm"]))
  expect_false(is.na(two_sided_mid_indices["Cpm"]))


  expect_lt(two_sided_off_indices[["Cpm"]], two_sided_mid_indices[["Cpm"]] + 1e-12)

  # Cpm cannot exceed Cp
  expect_lte(two_sided_off_indices[["Cpm"]], two_sided_off_indices[["Cp"]] + 1e-12)
  expect_lte(two_sided_mid_indices[["Cpm"]], two_sided_mid_indices[["Cp"]] + 1e-12)
})
