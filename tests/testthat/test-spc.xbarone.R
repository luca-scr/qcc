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

  testthat::describe("MSSD Estimator", {
    it("matches the analytic estimate for two observations", {
      # c4'(2) = sqrt(2 / pi), so sigma-hat = abs(diff(x)) * sqrt(pi) / 2.
      expect_equal(
        sd.xbar.one(c(1, 3), std.dev = "MSSD"),
        sqrt(pi),
        tolerance = 1e-12
      )
    })

    it("corrects successive squared differences using the MSSD factor", {
      sim <- readRDS(testthat::test_path("fixtures", "c4_mssd_mc.rds"))
      # Successive squared differences are 4, 1, 16, 4; n = 5.
      expect_equal(
        sd.xbar.one(c(1, 3, 2, 6, 4), std.dev = "MSSD"),
        sqrt(25 / 8) / sim$estimate[4],
        tolerance = 1e-6
      )
    })

    it("omits missing observations before forming successive differences", {
      expect_equal(
        sd.xbar.one(c(NA, 1, NA, NaN, 3, NA), std.dev = "MSSD"),
        sqrt(pi),
        tolerance = 1e-12
      )
    })

    it("returns NA when fewer than two observations remain", {
      for (x in list(numeric(), NA_real_, c(NA, NaN), 1, c(NA, 1, NA))) {
        expect_identical(sd.xbar.one(x, std.dev = "MSSD"), NA_real_)
      }
    })

    it("returns zero for constant observations", {
      expect_equal(sd.xbar.one(rep(4, 5), std.dev = "MSSD"), 0)
    })

    it("accepts a one-column matrix and ignores the MR window size", {
      expect_equal(
        sd.xbar.one(matrix(c(1, 3), ncol = 1), std.dev = "MSSD", r = 5),
        sqrt(pi),
        tolerance = 1e-12
      )
    })
  })
})
