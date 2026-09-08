data <- matrix(
  c(
    10.00, 10.02, 9.98, 10.01,
    10.03, 9.99, 10.01, 10.00,
    9.97, 10.04, 10.02, 9.99
  ),
  nrow = 3,
  byrow = TRUE
)

xbar <- ocCurves(
  qcc(data, type = "xbar", nsigmas = 3, plot = FALSE),
  size = 4,
  shift = c(0, 1, 2),
  nsigmas = 3
)
ranges <- ocCurves(
  qcc(data, type = "R", confidence.level = 0.95, plot = FALSE),
  size = 4,
  multiplier = c(1, 2, 3),
  nsigmas = NULL
)
standard_deviations <- ocCurves(
  qcc(data, type = "S", confidence.level = 0.95, plot = FALSE),
  size = 4,
  multiplier = c(1, 2, 3),
  nsigmas = NULL
)

test_that("continuous OC curves calculate non-detection probabilities", {
  expect_equal(
    unname(xbar$beta[, 1]),
    c(0.99730020, 0.84134446, 0.15865525),
    tolerance = 5e-7
  )

  expect_equal(
    unname(ranges$beta[, 1]),
    c(0.95000001, 0.50280230, 0.21521470),
    tolerance = 5e-7
  )

  expect_equal(
    unname(standard_deviations$beta[, 1]),
    c(0.95000000, 0.49127011, 0.20713356),
    tolerance = 5e-7
  )
})

test_that("continuous OC curves preserve fields and matrix labels", {
  for (type in c("xbar", "R", "S")) {
    grid_name <- if (type == "xbar") "shift" else "multiplier"
    axis_label <- if (type == "xbar") "shift (StdDev)" else "scale multiplier"
    chart <- qcc(data, type = type, confidence.level = 0.95, plot = FALSE)

    for (nsigmas in list(3, NULL)) {
      for (grid in list(1.25, c(1, 1.25, 2.5))) {
        for (size in list(4, c(8, 2, 4, 2))) {
          args <- list(object = chart, size = size, nsigmas = nsigmas)
          args[[grid_name]] <- grid
          curves <- do.call(ocCurves, args)
          expected_size <- if (length(size) == 1) 4 else c(2, 4, 8)
          expected_rows <- if (length(grid) == 1) "1.25" else c("1.00", "1.25", "2.50")
          expected_dimnames <- setNames(
            list(expected_rows, as.character(expected_size)),
            c(axis_label, "sample size")
          )

          expect_identical(class(curves), "ocCurves")
          expect_identical(names(curves), c("type", "size", grid_name, "beta", "ARL"))
          expect_identical(curves$type, type)
          expect_identical(curves$size, expected_size)
          expect_identical(curves[[grid_name]], grid)
          expect_identical(dimnames(curves$beta), expected_dimnames)
          expect_identical(curves$ARL, 1 / (1 - curves$beta))
        }
      }
    }
  }
})

test_that("binomial OC curves support p and np charts", {
  counts <- c(1, 2, 3, 4, 2)
  sizes <- rep(20, length(counts))

  expect_warning(
    proportions <- ocCurves(
      qcc(counts, sizes = sizes, type = "p", plot = FALSE)
    ),
    "discreteness of the binomial distribution"
  )
  expect_warning(
    counts_per_sample <- ocCurves(
      qcc(counts, sizes = sizes, type = "np", plot = FALSE)
    ),
    "discreteness of the binomial distribution"
  )

  expect_s3_class(proportions, "ocCurves")
  expect_s3_class(counts_per_sample, "ocCurves")
  for (curves in list(proportions, counts_per_sample)) {
    expect_identical(names(curves), c("type", "p", "beta", "ARL"))
    expected_rows <- sprintf("%.2f", seq(0, 1, length.out = 101))
    expect_identical(dimnames(curves$beta),
                     setNames(list(expected_rows, "beta"), c("fraction nonconforming", "")))
    expect_identical(dimnames(curves$ARL),
                     setNames(list(expected_rows, "ARL"), c("fraction nonconforming", "")))
  }
  expect_identical(proportions$type, "p")
  expect_identical(counts_per_sample$type, "np")
  expect_equal(proportions$p, seq(0, 1, length.out = 101))
  expect_equal(
    unname(proportions$beta[c(1, 6, 11, 51, 101), 1]),
    c(1, 0.99996605, 0.99761391, 0.05765915, 0),
    tolerance = 5e-8
  )
  expect_equal(counts_per_sample$beta, proportions$beta)
  expect_equal(
    unname(proportions$ARL),
    unname(1 / (1 - proportions$beta))
  )
  expect_equal(
    unname(counts_per_sample$ARL),
    unname(1 / (1 - counts_per_sample$beta))
  )
})

test_that("Poisson OC curves support c and u charts", {
  counts <- c(2, 3, 4, 3, 2)

  expect_warning(
    defects <- ocCurves(qcc(counts, type = "c", plot = FALSE)),
    "discreteness of the Poisson distribution"
  )
  expect_warning(
    defects_per_unit <- ocCurves(
      qcc(counts, sizes = rep(5, length(counts)), type = "u", plot = FALSE)
    ),
    "discreteness of the Poisson distribution"
  )

  expect_s3_class(defects, "ocCurves")
  expect_s3_class(defects_per_unit, "ocCurves")
  for (curves in list(defects, defects_per_unit)) {
    expect_identical(names(curves), c("type", "lambda", "beta", "ARL"))
    expect_identical(curves$lambda, 0:20)
    expected_rows <- paste0(0:20, ".00")
    expect_identical(dimnames(curves$beta),
                     setNames(list(expected_rows, "beta"), c("average nonconforming", "")))
    expect_identical(dimnames(curves$ARL),
                     setNames(list(expected_rows, "ARL"), c("average nonconforming", "")))
  }
  expect_identical(defects$type, "c")
  expect_identical(defects_per_unit$type, "u")
  expect_equal(
    unname(defects$beta[c(1, 4, 8, 13, 21), 1]),
    c(1, 0.98809550, 0.59871384, 0.08950450, 0.00077859),
    tolerance = 5e-8
  )
  expect_equal(defects_per_unit$beta, defects$beta)
  expect_equal(unname(defects$ARL), unname(1 / (1 - defects$beta)))
  expect_equal(
    unname(defects_per_unit$ARL),
    unname(1 / (1 - defects_per_unit$beta))
  )
})

test_that("OC curves reject unsupported chart types", {
  data <- matrix(c(10.1, 9.9, 10.0, 10.2, 9.8), ncol = 1)
  chart <- qcc(data, type = "xbar.one", plot = FALSE) # NOTE: Has to be supported by qcc

  expect_error(
    ocCurves(chart),
    "Operating characteristic curves not available for this type of chart"
  )
})
#
test_that("ocCurves matches its visual snapshot", {

  vdiffr::expect_doppelganger(
    "xbar OC curves",
    plot(xbar)
  )
  vdiffr::expect_doppelganger(
    "R OC curves",
    plot(ranges)
  )
  vdiffr::expect_doppelganger(
    "S OC curves",
    plot(standard_deviations)
  )
})
