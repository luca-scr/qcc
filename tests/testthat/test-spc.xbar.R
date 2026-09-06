# FIX: stats.xbar() and sd.xbar() do not expand a scalar `sizes` when called
#   directly, and mismatched lengths are recycled into incorrect estimates.
# FIX: stats.xbar() and sd.xbar() do not validate values in `sizes`. Noninteger
#   values are accepted, while nonpositive or missing values produce invalid
#   results or cryptic errors from the constants calculations.
# TODO: Define behavior for completely missing subgroups. stats.xbar() returns
#   NaN, while sd.xbar() fails after inferring a subgroup size of zero.
# TODO: Decide whether numeric `std.dev` must be scalar. sd.xbar() accepts a
#   vector and returns it unchanged, but the documentation describes one
#   process standard deviation.
# FIX: limits.xbar() does not validate `nsigmas`, `sizes`, or `std.dev`, so
#   invalid values can produce reversed, infinite, NaN, or recycled limits.
# FIX: limits.xbar() reports a base missing-value error for NA or NaN `conf`
#   instead of the documented invalid-conf error.
# PERFORMANCE: The MVLUE-R estimator calculates d2 for the same sizes three
#   times, including the calculation performed internally by d3.

testthat::describe("stats.xbar()", {
  data <- rbind(
    first = c(2,  4, 6,  8),
    second = c(1, 5, NA, NA),
    third = c(3,  6, 9,  NA)
  )
  sizes <- c(4, 2, 3)

  it("computes subgroup means and their size-weighted center", {
    result <- stats.xbar(data, sizes = sizes)

    expect_equal(result$statistics, c(first = 5, second = 3, third = 6))
    expect_equal(result$center, 44 / 9)
  })

  it("infers subgroup sizes from nonmissing observations", {
    expect_equal(
      stats.xbar(data),
      stats.xbar(data, sizes = sizes)
    )
  })
})


testthat::describe("sd.xbar()", {
  data <- rbind(
    c(2, 4,  6, 8),
    c(1, 5, NA, NA),
    c(3, 6,  9, NA)
  )
  sizes <- c(4, 2, 3)

  it("defaults to the UWAVE-R estimator", {
    expect_equal(
      sd.xbar(data),
      sd.xbar(data, sizes = sizes, std.dev = "UWAVE-R"),
      tolerance = testthat_tolerance()
    )
  })

  it("infers subgroup sizes from nonmissing observations", {
    methods <- c("UWAVE-R", "UWAVE-SD", "MVLUE-R", "MVLUE-SD", "RMSDF")
    inferred <- vapply(methods, function(method) {
      sd.xbar(data, std.dev = method)
    }, numeric(1))
    supplied <- vapply(methods, function(method) {
      sd.xbar(data, sizes = sizes, std.dev = method)
    }, numeric(1))

    expect_equal(inferred, supplied, tolerance = testthat_tolerance())
  })

  testthat::describe("range estimators", {
    it("reproduces UWAVE-R estimates", {
      expect_equal(
        sd.xbar(data, std.dev = "UWAVE-R"),
        3.334735,
        tolerance = 5e-4
      )
    })

    it("reproduces MVLUE-R estimates", {
      expect_equal(
        sd.xbar(data, std.dev = "MVLUE-R"),
        3.226929,
        tolerance = 5e-4
      )
    })
  })

  testthat::describe("standard-deviation estimators", {
    it("reproduces UWAVE-SD estimates", {
      expect_equal(
        sd.xbar(data, std.dev = "UWAVE-SD"),
        3.244180,
        tolerance = 5e-4
      )
    })

    it("reproduces MVLUE-SD estimates", {
      expect_equal(
        sd.xbar(data, std.dev = "MVLUE-SD"),
        3.113833,
        tolerance = 5e-4
      )
    })

    it("reproduces RMSDF estimates", {
      expect_equal(
        sd.xbar(data, std.dev = "RMSDF"),
        2.886142,
        tolerance = 5e-4
      )
    })
  })

  it("returns a supplied numeric standard deviation", {
    expect_identical(
      sd.xbar(data, sizes = sizes, std.dev = 2.5),
      2.5
    )
  })

  # TODO: Check correctness; some estimators are weighted
  it("rejects subgroups with one observation", {
    expect_error(
      sd.xbar(rbind(c(1, NA), c(2, 3))),
      "group sizes must be larger than one"
    )
  })

  it("rejects unknown estimators", {
    expect_error(
      sd.xbar(data, std.dev = "unknown"),
      "should be one of"
    )
  })
})


testthat::describe("limits.xbar()", {
  center <- 10
  std.dev <- 2
  sizes <- c(2, 4, 8)

  it("requires either sigma or probability limits", {
    expect_error(
      limits.xbar(center, std.dev, sizes),
      "Argument 'nsigmas' or 'conf' must be provided"
    )
  })

  testthat::describe("sigma limits", {

    it("computes limits for each distinct subgroup size", {
      expect_equal(
        unname(limits.xbar(center, std.dev, sizes, nsigmas = 3)),
        cbind(
          c(5.757359, 7.000000, 7.878680),
          c(14.24264, 13.00000, 12.12132)
        ),
        tolerance = 5e-4
      )
    })

    it("returns one limits row for a common subgroup size", {
      limits <- limits.xbar(
        center,
        std.dev,
        sizes = rep(4, 3),
        nsigmas = 3
      )

      expect_equal(dim(limits), c(1L, 2L))
      expect_equal(unname(limits), matrix(c(7, 13), nrow = 1))
    })
  })

  testthat::describe("probability limits", {

    # HACK: impossible to test without reimplementation.
    it("uses the two-sided normal quantile", {
      conf <- 0.95

      expect_equal(
        unname(limits.xbar(center, std.dev, sizes, conf = 0.95)),
        cbind(
          c(7.228192, 8.040036, 8.614096),
          c(12.77181, 11.95996, 11.38590)
        ),
        tolerance = 5e-4
      )
    })

    it("rejects invalid confidence levels", {
      invalid <- list(-0.1, 0, 1, 1.1, "0.95", c(0.9, 0.95))

      for (conf in invalid) {
        expect_error(
          limits.xbar(center, std.dev, sizes, conf = conf),
          "invalid 'conf' argument"
        )
      }
    })
  })
})
