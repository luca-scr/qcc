testthat::describe("mssd_bias_correction_factor", {
  it("returns tabulated values", {
    # reference table from Minitab: https://support.minitab.com/en-us/minitab/help-and-how-to/quality-and-process-improvement/capability-analysis/how-to/capability-analysis/normal-capability-analysis/methods-and-formulas/methods/#unbiasing-constant-c4
    tabulated_values <- readRDS("tests/testthat/fixtures/mssd_bias_correction.rds")
    expect_equal(.c4p(2:500), tabulated_values[-1], tolerance = 5e-7)
  })
})

