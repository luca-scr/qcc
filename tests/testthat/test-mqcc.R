make_mqcc_phase1_data <- function() {
  x1 <- matrix(
    c(
      10, 11, 9,
      11, 10, 12,
      9, 10, 11,
      10, 9, 12
    ),
    nrow = 4,
    byrow = TRUE
  )
  x2 <- matrix(
    c(
      20, 19, 21,
      18, 20, 19,
      22, 21, 23,
      19, 18, 20
    ),
    nrow = 4,
    byrow = TRUE
  )

  list(X1 = x1, X2 = x2)
}

make_mqcc_newdata <- function() {
  x1 <- matrix(
    c(
      10, 10, 11,
      12, 11, 10
    ),
    nrow = 2,
    byrow = TRUE
  )
  x2 <- matrix(
    c(
      20, 21, 19,
      18, 19, 20
    ),
    nrow = 2,
    byrow = TRUE
  )

  list(X1 = x1, X2 = x2)
}

test_that("mqcc constructor returns class, type, and phase-I labels", {
  phase1 <- make_mqcc_phase1_data()

  chart <- mqcc(
    phase1,
    plot = FALSE,
    limits = FALSE,
    pred.limits = FALSE
  )

  expect_s3_class(chart, "mqcc")
  expect_equal(chart$type, "T2")
  expect_equal(chart$var.names, c("X1", "X2"))
  expect_equal(names(chart$statistics), as.character(seq_len(nrow(phase1$X1))))
})

test_that("mqcc newdata stores newstats labels and preserves phase-I labels", {
  phase1 <- make_mqcc_phase1_data()
  phase2 <- make_mqcc_newdata()
  phase1_labels <- c("phase-1", "phase-2", "phase-3", "phase-4")

  with_supplied_labels <- mqcc(
    phase1,
    labels = phase1_labels,
    newdata = phase2,
    newlabels = c("new-1", "new-2"),
    plot = FALSE,
    limits = FALSE,
    pred.limits = FALSE
  )

  expect_equal(names(with_supplied_labels$statistics), phase1_labels)
  expect_equal(names(with_supplied_labels$newstats), c("new-1", "new-2"))

  with_generated_labels <- mqcc(
    phase1,
    labels = phase1_labels,
    newdata = phase2,
    plot = FALSE,
    limits = FALSE,
    pred.limits = FALSE
  )

  expected_generated <- as.character(
    seq.int(
      length(with_generated_labels$statistics) + 1,
      length(with_generated_labels$statistics) + length(with_generated_labels$newstats)
    )
  )

  expect_equal(names(with_generated_labels$statistics), phase1_labels)
  expect_equal(names(with_generated_labels$newstats), expected_generated)

})

test_that("mqcc validates confidence.level and newlabels length", {
  phase1 <- make_mqcc_phase1_data()
  phase2 <- make_mqcc_newdata()

  expect_error(
    mqcc(phase1, confidence.level = 1, plot = FALSE),
    "confidence.level must be a numeric value in the range \\(0,1\\)"
  )

  expect_error(
    mqcc(
      phase1,
      newdata = phase2,
      newlabels = "new-1",
      plot = FALSE
    ),
    "labels must match the length of samples provided"
  )
})

test_that("no visual regressions in ellipseChart", {
  phase1 <- make_mqcc_phase1_data()
  phase2 <- make_mqcc_newdata()

  chart <- mqcc(
    phase1,
    labels = c("old-1", "old-2", "old-3", "old-4"),
    newdata = phase2,
    newlabels = c("new-1", "new-2"),
    plot = FALSE,
  )

  vdiffr::expect_doppelganger(
    "mqcc plot",
    plot(chart)
  )

  vdiffr::expect_doppelganger(
    "mqcc ellipseChart",
    ellipseChart(chart)
  )
})
