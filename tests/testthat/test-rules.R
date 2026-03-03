make_rule_test_qcc <- function(statistics, newstats = numeric(), center = 0, std.dev = 1, nsigmas = 3) {
  limits <- matrix(c(center - nsigmas * std.dev, center + nsigmas * std.dev), ncol = 2)
  colnames(limits) <- c("LCL", "UCL")

  structure(
    list(
      statistics = statistics,
      newstats = newstats,
      limits = limits,
      center = center,
      std.dev = std.dev,
      nsigmas = nsigmas,
      type = "xbar.one",
      sizes = rep(1, length(statistics)),
      newsizes = rep(1, length(newstats))
    ),
    class = "qcc"
  )
}

test_that("WER1-WER4 return expected violating indices for deterministic qcc objects", {
  wer1_object <- make_rule_test_qcc(c(0, 3.2, -3.5, 2.9))
  expect_equal(qccRulesViolatingWER1(wer1_object), c(2, 3))

  wer2_object <- make_rule_test_qcc(c(0, 2.5, 1.0, 2.6, 0))
  expect_equal(qccRulesViolatingWER2(wer2_object), 4)

  wer3_object <- make_rule_test_qcc(c(0, 1.2, 1.3, 0.2, 1.4, 1.5))
  expect_equal(qccRulesViolatingWER3(wer3_object), 6)

  wer4_object <- make_rule_test_qcc(c(-1, -1, -1, -1, -1, -1, -1, -1, 1))
  expect_equal(qccRulesViolatingWER4(wer4_object), 8)
})

test_that("qccRules precedence reports WER1 when overlaps include multiple rules", {
  object <- make_rule_test_qcc(c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 3.2))

  expect_equal(qccRulesViolatingWER4(object), 8)
  expect_equal(qccRulesViolatingWER3(object), 5:8)
  expect_equal(qccRulesViolatingWER2(object), 3:8)
  expect_equal(qccRulesViolatingWER1(object), 8)

  wer4_only <- qccRules(object, rules = 4)
  expect_equal(wer4_only[8], 4)
  expect_true(all(is.na(wer4_only[1:7])))

  result_432 <- qccRules(object, rules = c(4, 3, 2))
  expect_equal(result_432, c(NA, NA, 2, 2, 2, 2, 2, 2), ignore_attr = TRUE)
  expect_equal(result_432[8], 2)

  result_4321 <- qccRules(object, rules = c(4, 3, 2, 1))
  expect_equal(result_4321, c(NA, NA, 2, 2, 2, 2, 2, 1), ignore_attr = TRUE)
  expect_equal(result_4321[8], 1)
})

test_that("qccRules errors for objects that are not qcc or mqcc", {
  expect_error(
    qccRules(list(statistics = 1), rules = 1),
    "input object must be of class 'qcc' or 'mqcc'"
  )
})
