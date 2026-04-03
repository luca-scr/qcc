INDICES <- c("Cp", "Cp_l", "Cp_u", "Cp_k", "Cpm", "Pp", "Pp_l", "Pp_u", "Pp_k", "Ppm")
INDEX_COLUMNS <- c("Value", "2.5%", "97.5%")

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
  expect_true(is.matrix(capability$indices))
  expect_equal(rownames(capability$indices), INDICES)
  expect_equal(colnames(capability$indices), INDEX_COLUMNS)
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
      1.372667, 0.808460, 1.937713,
      1.441300, 0.911657, 1.970943,
      1.304033, 0.820114, 1.787953,
      1.304033, 0.727408, 1.880659,
      1.344463, 0.804207, 1.885321,
      1.612598, 0.949773, 2.276410,
      1.693228, 1.078708, 2.307749,
      1.531968, 0.971902, 2.092035,
      1.531968, 0.864608, 2.199329,
      1.567395, 0.933327, 2.202237
    ),
    nrow = length(INDICES),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(
      INDICES,
      INDEX_COLUMNS
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
      NA_real_, NA_real_, NA_real_,
      NA_real_, NA_real_, NA_real_,
      1.304033, 0.820114, 1.787953,
      1.304033, 0.727408, 1.880659,
      NA_real_, NA_real_, NA_real_,
      NA_real_, NA_real_, NA_real_,
      NA_real_, NA_real_, NA_real_,
      1.531968, 0.971902, 2.092035,
      1.531968, 0.864608, 2.199329,
      NA_real_, NA_real_, NA_real_
    ),
    nrow = length(INDICES),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(
      INDICES,
      INDEX_COLUMNS
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
      NA_real_, NA_real_, NA_real_,
      1.441300, 0.911657, 1.970943,
      NA_real_, NA_real_, NA_real_,
      1.441300, 0.810191, 2.072409,
      NA_real_, NA_real_, NA_real_,
      NA_real_, NA_real_, NA_real_,
      1.693228, 1.078708, 2.307749,
      NA_real_, NA_real_, NA_real_,
      1.693228, 0.960982, 2.425474,
      NA_real_, NA_real_, NA_real_
    ),
    nrow = length(INDICES),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(
      INDICES,
      INDEX_COLUMNS
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
      1.372667, 0.808460, 1.937713,
      1.441300, 0.911657, 1.970943,
      1.304033, 0.820114, 1.787953,
      1.304033, 0.727408, 1.880659,
      1.167834, 0.650542, 1.686689,
      1.612598, 0.949773, 2.276410,
      1.693228, 1.078708, 2.307749,
      1.531968, 0.971902, 2.092035,
      1.531968, 0.864608, 2.199329,
      1.305161, 0.712493, 1.899988
    ),
    nrow = length(INDICES),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(
      INDICES,
      INDEX_COLUMNS
    )
  )

  expect_equal(
    two_sided_off_capability$indices,
    expected_indices,
    tolerance = 1e-6
  )
})
