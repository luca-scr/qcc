# TODO: Add describe blocks for stats.xbar.one and limits.xbar.one

# TODO: Specify what it does with imputs that are not atomic vectors
# TODO: argument `r` is not validated or tested
# FIX: Ignoring NAs should be handled at the sd.xbar.one level, not implemented
#   for each estimator. But how would it handle Estimators of autocorrelated
#   individual observations?
testthat::describe("sd.xbar.one", {

  it("defaults to the MR estimator on successive pairs",{
    expect_equal(
      sd.xbar.one(antifreeze$water),
      sd.xbar.one(antifreeze$water, std.dev = "MR", r = 2),
      tolerance = testthat_tolerance()
    )
  })

  testthat::describe("MR Estimator",{
    it("can reproduce estimates", {
      expect_equal(
        sd.xbar.one(antifreeze$water, std.dev = "MR"),
        0.179393,
        tolerance = 5e-4
      )
    })

    # This behavior should be optional, because it is only correct under the
    #   assumption that missing data and non-missing data share the same
    #   distribution. or warn when omitting NAs
    it("ignores missing data", {
      mat <- rbind(
        c(100, 110, 102, 99, NA,  105, 98, 112),
        c(100, 110, 102, 99, NA,   NA, 98, 112),
        c(NA,  110, 102, 99, NA,  105, 98, NA),
        c(NA,   NA, 102, 99, NA,  105, 98, 112)
      )
      expect_equal(
        apply(mat, 1, sd.xbar.one),
        c(7.089815, 6.380833, 5.317361, 6.646701),
        tolerance = 5e-4
      )
    })
  })

  testthat::describe("SD Estimator", {
    it("can reproduce estimates", {
      expect_equal(
        sd.xbar.one(antifreeze$water, std.dev = "SD"),
        0.221679,
        tolerance = 5e-4
      )
    })
    it("ignores missing data", {
      mat <- rbind(
        c(100, 110, 102, 99, NA,  105, 98, 112),
        c(100, 110, 102, 99, NA,   NA, 98, 112),
        c(NA,  110, 102, 99, NA,  105, 98, NA),
        c(NA,   NA, 102, 99, NA,  105, 98, 112)
      )
      expect_equal(
        apply(mat, 1, sd.xbar.one, std.dev = "SD"),
        c(5.731807, 6.296851, 5.179084, 5.989746),
        tolerance = 5e-4
      )
    })
  })
})
