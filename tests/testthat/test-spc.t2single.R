# TODO: Define missing-value behavior. The default center and covariance
#   propagate missing values to every statistic without an informative error.
# FIX: limits.T2.single() does not validate confidence levels or distribution
#   degrees of freedom, so invalid inputs can produce NaN or infinite limits.
# TODO: Validate that `size` is one or remove it; limits.T2.single() currently
#   accepts and ignores every supplied value.
# FIX: The @return documentation calls the covariance component `S`, but the
#   returned component is named `cov`.
# FIX: the following tests assert a stable implementation. Not tests for
#   correctness!

testthat::describe("stats.T2.single()", {
  it("reproduces individual-observation Hotelling statistics", {
    expect_snapshot(stats.T2.single(boiler))
  })

  # TODO: test user-supplied center and cov
})

testthat::describe("limits.T2.single()", {
  it("reproduces phase I and phase II limits", {
    expect_snapshot(
      limits.T2.single(ngroups = 25, nvars = 8, conf = 0.999)
    )
  })
})
