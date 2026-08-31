# HACK: We use `testthat::describe` because we have `qcc::describe`
testthat::describe("d2()", {
  it("returns tabulated values", {
    expected <- c( # SOURCE: d2 <- qcc::qcc.options("exp.R.unscaled")
      1.128, 1.693, 2.059, 2.326, 2.534, 2.704, 2.847, 2.970,
      3.078, 3.173, 3.258, 3.336, 3.407, 3.472, 3.532, 3.588,
      3.640, 3.689, 3.735, 3.778, 3.819, 3.858, 3.895, 3.931
    )
    expect_equal(d2(2:25), expected, tolerance = 5e-4)
  })

  it("reasonably approximates Wardell2025 analytic solutions", {
    expect_equal(
      d2(2:5),
      c(
        2 / sqrt(pi),
        3 / sqrt(pi),
        12 / (pi * sqrt(pi)) * atan(sqrt(2)),
        30 / pi^(3 / 2) * atan(sqrt(2)) - 5 / sqrt(pi)
      ),
      tolerance = testthat_tolerance()
    )
  })

  test_that("Handles unsupported sample sizes", {
    expect_warning(x <- d2(c(0, 1, 20.5, Inf, -Inf, NA, NaN, 2)))
    expect_equal(x, c(rep(NA, 7), d2(2)))
  })
})


testthat::describe("d3()", {
  it("returns tabulated values",{
    expected <- c( # SOURCE: d3 <- qcc::qcc.options("se.R.unscaled")
      0.8525033, 0.8883697, 0.8798108, 0.8640855, 0.8480442, 0.8332108, 0.8198378,
      0.8078413, 0.7970584, 0.7873230, 0.7784873, 0.7704257, 0.7630330, 0.7562217,
      0.7499188, 0.7440627, 0.7386021, 0.7334929, 0.7286980, 0.7241851, 0.7199267,
      0.7158987, 0.7120802, 0.7084528, 0.7050004, 0.7017086, 0.6985648, 0.6955576,
      0.6926770, 0.6899137, 0.6872596, 0.6847074, 0.6822502, 0.6798821, 0.6775973,
      0.6753910, 0.6732584, 0.6711952, 0.6691976, 0.6672619, 0.6653848, 0.6635632,
      0.6617943, 0.6600754, 0.6584041, 0.6567780, 0.6551950, 0.6536532, 0.6521506
    )
    expect_equal(d3(2:50), expected, tolerance = 5e-4)

    # Test input vectors with repeated elements because our implementation [integrate_ok()]
    #   deduplicates its input then maps back to the original vector length and order,
    #   for the sake of better performance.
    expect_equal(
      d3(c(rep(2,7),5, rep(6,7))), 
      expected[c(rep(1, 7), 4, rep(5, 7))],
      tolerance = 5e-4
    )
  })

  it("reasonably approximates Wardell2025 analytic solutions", {
    expect_equal(
      d3(2:5),
      c(
        sqrt(2 - 4 / pi),
        sqrt((2 * pi + 3 * sqrt(3) - 9) / pi),
        sqrt((2 * pi + 2 * sqrt(3) + 6) / pi - d2(4)^2),
        sqrt((2 * pi^2 + 10 * sqrt(3) * (atan(sqrt(5 / 3)) + 2 * sqrt(3) * atan(sqrt(1 / 5)))) / pi^2 - d2(5)^2)
      ),
      tolerance = testthat_tolerance()
    )
  })

  it("handles unsupported sample sizes", {
    expect_warning(x <- d3(c(0, 1, 20.5, Inf, -Inf, NA, NaN, 2)))
    expect_equal(x, c(rep(NA, 7), d3(2)))
  })
})

testthat::describe("c4()", {
  it("returns expected values",{
    expected <- c( # SOURCE: c4 <- qcc:::qcc.c4(2:50)
      0.7978846, 0.8862269, 0.9213177, 0.9399856, 0.9515329, 0.9593688, 0.9650305,
      0.9693107, 0.9726593, 0.9753501, 0.9775594, 0.9794056, 0.9809714, 0.9823162,
      0.9834835, 0.9845064, 0.9854100, 0.9862141, 0.9869343, 0.9875829, 0.9881703,
      0.9887045, 0.9891927, 0.9896404, 0.9900525, 0.9904330, 0.9907856, 0.9911130,
      0.9914181, 0.9917028, 0.9919693, 0.9922192, 0.9924540, 0.9926751, 0.9928836,
      0.9930805, 0.9932668, 0.9934434, 0.9936109, 0.9937701, 0.9939216, 0.9940659,
      0.9942034, 0.9943348, 0.9944603, 0.9945804, 0.9946954, 0.9948056, 0.9949113
    )
    expect_equal(c4(2:50), expected, tolerance = 5e-4)
  })

  it("handles unsupported sample sizes", {
    expect_warning(x <- c4(c(0, 1, 20.5, Inf, -Inf, NA, NaN, 2)))
    expect_equal(x, c(rep(NA, 7), c4(2)))
  })
})
