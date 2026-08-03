# FIX: stats.T2() silently recycles unequal subgroup sizes when called directly.
#   mqcc() rejects them, but stats.T2() is exported and performs no validation.
# TODO: Define missing-value behavior. Means and centers ignore missing values,
#   while the default pooled covariance propagates them to every statistic.
# FIX: limits.T2() does not validate confidence levels or distribution degrees
#   of freedom, so invalid inputs can produce NaN or infinite limits.
# FIX: The @return documentation calls the covariance component `S`, but the
#   returned component is named `cov`.
# FIX: the following tests assert a stable implementation. Not tests for
#   correctness!
# HACK: we use expect_snapshot instead of expect_equal because expected output is large.

testthat::describe("stats.T2()", {
  it("reproduces grouped Hotelling statistics", {
    expect_snapshot(stats.T2(RyanMultivar))
  })

  it("uses supplied center and covariance estimates", {
    center <- c(X1 = 65, X2 = 20)
    cov <- matrix(
      c(200, 105, 105, 59),
      nrow = 2,
      dimnames = list(c("X1", "X2"), c("X1", "X2"))
    )

    expect_snapshot(
      stats.T2(RyanMultivar, center = center, cov = cov)
    )
  })
})


testthat::describe("limits.T2()", {
  it("reproduces phase I and phase II limits", {
    expect_snapshot(
      limits.T2(
        ngroups = 20,
        size = 4,
        nvars = 2,
        conf = (1 - 0.0027)^2
      )
    )
  })
})
