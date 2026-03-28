INDICES <- c("Cp", "Cp_l", "Cp_u", "Cp_k", "Cpm")

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
two_sided_specs <- c(9.9, 10.1)
midpoint_target <- mean(two_sided_specs)
off_target <- 10.02
chart <- qcc(make_grouped_xbar_data(), type = "xbar", nsigmas = 3)

upper_capability <- processCapability(
  chart,
  spec.limits = c(NA_real_, two_sided_specs[2]),
  target = off_target
)
lower_capability <- processCapability(
  chart,
  spec.limits = c(two_sided_specs[1], NA_real_),
  target = off_target
)
two_sided_off_capability <- processCapability(
  chart,
  spec.limits = two_sided_specs,
  target = off_target
)
two_sided_mid_capability <- processCapability(
  chart,
  spec.limits = two_sided_specs,
  target = midpoint_target
)

capability <- two_sided_mid_capability

test_that("ProcessCapability returns expected structure", {
  expect_s3_class(capability, "processCapability")
  expect_equal(names(capability$spec.limits), c("LSL", "USL"))
  expect_equal(unname(capability$spec.limits), two_sided_specs)
})

test_that("processCapability returns indices matrix with expected structure", {
  expect_true(is.matrix(capability$indices))
  expect_equal(rownames(capability$indices), INDICES)
  expect_equal(colnames(capability$indices)[1], "Value")
  expect_equal(ncol(capability$indices), 3)
})

test_that("ProcessCapability rejects invalid specs", {
  expect_error(
    processCapability(chart, spec.limits = c(NA_real_, NaN)),
    "invalid specification limits"
  )
})

test_that("ProcessCapability warns and defaults target to the midpoint when target is omitted", {
  expect_message(
    capability <- processCapability(chart, spec.limits = two_sided_specs),
    "target value not provided; using midpoint of specification limits"
  )
  expect_true(capability$has.target)
  expect_equal(capability$target, mean(two_sided_specs))
})

test_that("ProcessCapability returns expected indices. Two-sided on-target case", {
  expected_indices <- matrix(
    c(
      1.372667, 1.441300, 1.304033, 1.304033, 1.344463,
      0.808460, 0.911657, 0.820114, 0.727408, 0.804207,
      1.937713, 1.970943, 1.787953, 1.880659, 1.885321
    ),
    nrow = 5,
    ncol = 3,
    dimnames = list(
      INDICES,
      c("Value", "2.5%", "97.5%")
    )
  )

  expect_equal(
    two_sided_mid_capability$indices,
    expected_indices,
    tolerance = 1e-6
  )
})

test_that("ProcessCapability returns expected indices. Upper spec only case", {
  expected_upper_indices <- matrix(
    c(
      NA_real_, NA_real_, 1.304033, 1.304033, NA_real_,
      NA_real_, NA_real_, 0.820114, 0.727408, NA_real_,
      NA_real_, NA_real_, 1.787953, 1.880659, NA_real_
    ),
    nrow = 5,
    ncol = 3,
    dimnames = list(
      INDICES,
      c("Value", "2.5%", "97.5%")
    )
  )

  expect_equal(
    upper_capability$indices,
    expected_upper_indices,
    tolerance = 1e-6
  )
})

test_that("ProcessCapability returns expected indices. Lower spec only case", {
  expected_lower_indices <- matrix(
    c(
      NA_real_, 1.441300, NA_real_, 1.441300, NA_real_,
      NA_real_, 0.911657, NA_real_, 0.810191, NA_real_,
      NA_real_, 1.970943, NA_real_, 2.072409, NA_real_
    ),
    nrow = 5,
    ncol = 3,
    dimnames = list(
      INDICES,
      c("Value", "2.5%", "97.5%")
    )
  )
  
  expect_equal(
    lower_capability$indices,
    expected_lower_indices,
    tolerance = 1e-6
  )
})

test_that("ProcessCapability returns expected indices. Two-sided off-target case", {
  expected_indices <- matrix(
    c(
      1.372667, 1.441300, 1.304033, 1.304033, 1.167834,
      0.808460, 0.911657, 0.820114, 0.727408, 0.650542,
      1.937713, 1.970943, 1.787953, 1.880659, 1.686689
    ),
    nrow = 5,
    ncol = 3,
    dimnames = list(
      INDICES,
      c("Value", "2.5%", "97.5%")
    )
  )

  expect_equal(
    two_sided_off_capability$indices,
    expected_indices,
    tolerance = 1e-6
  )
})

upper_indices <- upper_capability$indices[, "Value"]
lower_indices <- lower_capability$indices[, "Value"]
two_sided_off_indices <- two_sided_off_capability$indices[, "Value"]
two_sided_mid_indices <- two_sided_mid_capability$indices[, "Value"]

test_that("Cp, Cp_l, Cp_u, Cp_k do not depend on target", {
  cp_family <- c("Cp", "Cp_l", "Cp_u", "Cp_k")

  expect_false(anyNA(two_sided_off_indices[cp_family]))
  expect_false(anyNA(two_sided_mid_indices[cp_family]))
  expect_equal(
    unname(two_sided_off_indices[cp_family]),
    unname(two_sided_mid_indices[cp_family]),
    tolerance = 1e-6
  )
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
})

test_that("Cpk collapses correctly to Cpu or Cpl in one-sided specs", {
  expect_equal(upper_indices[["Cp_k"]], upper_indices[["Cp_u"]], tolerance = 1e-6)
  expect_equal(lower_indices[["Cp_k"]], lower_indices[["Cp_l"]], tolerance = 1e-6)
})

test_that("Cpm calculation passes mathematical bounds", {
  expect_lt(two_sided_off_indices[["Cpm"]], two_sided_mid_indices[["Cpm"]] + 1e-12)
  expect_lte(two_sided_off_indices[["Cpm"]], two_sided_off_indices[["Cp"]] + 1e-12)
  expect_lte(two_sided_mid_indices[["Cpm"]], two_sided_mid_indices[["Cp"]] + 1e-12)
})
