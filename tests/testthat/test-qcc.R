make_grouped <- function() {
  matrix(
    c(
      10, 11, 12, 13,
      11, 10, 12, 14,
      9, 10, 11, 13
    ),
    nrow = 3,
    byrow = TRUE
  )
}

make_u_chart <- function() {
  data(circuit, package = "qcc")
  outofctrl <- c(6, 20)

  with(
    circuit[-outofctrl, ],
    qcc::qcc(
      x[trial],
      sizes = size[trial],
      type = "u",
      newdata = x[!trial],
      newsizes = size[!trial],
      rules = c(1, 2, 3)
    )
  )
}

test_that("Creating xbar charts return qcc object with two control limits", {
  data("pistonrings", package = "qcc")
  diameter <- qcc::qccGroups(data = pistonrings, diameter, sample)

  chart <- qcc::qcc(diameter[1:25, ], type = "xbar")

  expect_s3_class(chart, "qcc")
  expect_equal(chart$type, "xbar")
  expect_equal(ncol(chart$limits), 2)
})

test_that("xbar.one charts produce the correct number of violations", {
  data(viscosity, package = "qcc")
  chart <- with(viscosity, qcc::qcc(viscosity[trial], type = "xbar.one"))

  expect_equal(sum(!is.na(chart$violations)), 1)
})

test_that("u charts with newdata detect expected violations", {
  data(circuit, package = "qcc")
  outofctrl <- c(6, 20)

  chart <- with(
    circuit[-outofctrl, ],
    qcc(
      x[trial],
      sizes = size[trial],
      type = "u",
      newdata = x[!trial],
      newsizes = size[!trial]
    )
  )

  expect_s3_class(chart, "qcc")
  expect_equal(chart$type, "u")
  expect_equal(ncol(chart$limits), 2)
  expect_equal(sum(!is.na(chart$violations)), 1)
})

test_that("xbar end-to-end covers print, summary, plot, and limits-mode switch", {
  LINES_PRODUCED_BY_PRINT = 12
  LINES_PRODUCED_BY_SUMMARY = 12

  data(pistonrings, package = "qcc")
  diameter <- qcc::qccGroups(data = pistonrings, diameter, sample)

  chart <- qcc::qcc(
    data = diameter[1:25, ],
    type = "xbar",
    confidence.level = 0.95
  )

  printed <- capture.output(print(chart))
  summarized <- capture.output(summary(chart))
  plotted <- qcc::plot.qcc(chart, add.stats = FALSE, fill = FALSE)

  expect_length(printed, LINES_PRODUCED_BY_PRINT)
  expect_length(summarized, LINES_PRODUCED_BY_SUMMARY)
  expect_equal(nrow(plotted$data), length(chart$statistics))
  expect_no_error(ggplot2::ggplot_build(plotted))

  expect_warning(
    chart_with_limits <- qcc::qcc(
      data = diameter[1:25, ],
      type = "xbar",
      std.dev = 1,
      limits = c(73.98, 74.02)
    ),
    "std.dev' is not used when limits is given"
  )
  expect_equal(unname(as.numeric(chart_with_limits$limits[1, ])), c(73.98, 74.02))
})

test_that("plot.qcc chart.all controls included groups and builds cleanly", {
  chart <- make_u_chart()
  total_groups <- length(chart$statistics) + length(chart$newstats)

  plot_all <- qcc::plot.qcc(
    chart,
    chart.all = TRUE,
    fill = TRUE,
    add.stats = FALSE
  )
  plot_new <- qcc::plot.qcc(
    chart,
    chart.all = FALSE,
    fill = TRUE,
    add.stats = FALSE
  )

  expect_equal(nrow(plot_all$data), total_groups)
  expect_equal(nrow(plot_new$data), length(chart$newstats))
  expect_no_error(ggplot2::ggplot_build(plot_all))
  expect_no_error(ggplot2::ggplot_build(plot_new))
})

test_that("plot.qcc accepts Date and POSIXct xtime values", {
  chart <- make_u_chart()
  total_groups <- length(chart$statistics) + length(chart$newstats)

  date_index <- as.Date("2024-01-01") + seq_len(total_groups) - 1
  expect_no_error(
    plot_with_stats <- qcc::plot.qcc(
      chart,
      xtime = date_index,
      chart.all = TRUE,
      fill = FALSE,
      add.stats = TRUE
    )
  )
  expect_true(any(class(plot_with_stats) %in% c("patchwork", "ggplot", "gg")))

  datetime_index <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC") +
    3600 * seq_len(total_groups)
  expect_no_error(
    qcc::plot.qcc(
      chart,
      xtime = datetime_index,
      chart.all = FALSE,
      fill = FALSE,
      add.stats = FALSE
    )
  )
})

test_that("plot.qcc validates xtime class", {
  chart <- make_u_chart()
  total_groups <- length(chart$statistics) + length(chart$newstats)

  expect_error(
    qcc::plot.qcc(chart, xtime = rep("bad", total_groups)),
    "xtime must be"
  )
})


test_that("S-chart helpers keep expected structure and non-negative lower limits", {
  grouped <- make_grouped()
  subgroup_sizes <- rep(4, nrow(grouped))

  s_stats <- qcc::stats.S(grouped, subgroup_sizes)
  s_sd <- qcc::sd.S(grouped, subgroup_sizes)
  s_limits <- qcc::limits.S(
    center = s_stats$center,
    std.dev = s_sd,
    sizes = subgroup_sizes,
    nsigmas = 3
  )

  expect_equal(length(s_stats$statistics), nrow(grouped))
  expect_length(s_limits, 2)
  expect_true(all(s_limits[, "LCL"] >= 0))
})

test_that("qcc validates required chart inputs and dimensions", {
  grouped <- make_grouped()

  expect_error(qcc::qcc(c(1, 2, 3), type = "u"), "sample 'sizes' must be given")
  expect_error(qcc::qcc(grouped, type = "not_a_chart"), "invalid")
  expect_error(
    qcc::qcc(grouped, type = "xbar", std.dev = list(1)),
    "std.dev"
  )
  expect_error(
    qcc::qcc(grouped, type = "xbar", sizes = c(4, 4)),
    "sizes length doesn't match"
  )
  expect_error(
    qcc::qcc(grouped, type = "xbar", limits = "bad"),
    "'limits' must be a vector of length 2 or a 2-columns matrix"
  )
  expect_error(
    qcc::qcc(
      grouped,
      type = "xbar",
      newdata = grouped[1:2, ],
      newsizes = c(4, 4, 4)
    ),
    "newsizes length doesn't match with newdata"
  )
  expect_error(
    qcc::qcc(
      c(1, 2, 3),
      type = "u",
      sizes = c(10, 10, 10),
      newdata = c(1, 2)
    ),
    "newsizes"
  )
})

test_that("c-chart helpers enforce unit sizes and compute sd from mean count", {
  expect_error(qcc::stats.c(c(1, 2, 3), sizes = c(1, 2, 1)), "all sizes")
  expect_equal(qcc::sd.c(c(1, 2, 3), sizes = rep(1, 3)), sqrt(mean(c(1, 2, 3))))
})

test_that("p and np helpers remain internally consistent", {
  counts <- c(2, 4, 3, 5)
  sample_sizes <- c(10, 12, 11, 13)

  p_stats <- qcc::stats.p(counts, sample_sizes)
  np_stats <- qcc::stats.np(counts, sample_sizes)
  p_sd <- qcc::sd.p(counts, sample_sizes)
  np_sd <- qcc::sd.np(counts, sample_sizes)
  p_limits <- qcc::limits.p(p_stats$center, p_sd, sample_sizes, nsigmas = 3)
  np_limits <- qcc::limits.np(np_stats$center, np_sd, sample_sizes, nsigmas = 3)

  expect_equal(np_stats$statistics / sample_sizes, p_stats$statistics)
  expect_equal(np_limits / sample_sizes, p_limits, tolerance = 1e-12)
})

test_that("R-chart helpers keep expected structure and enforce max subgroup size", {
  grouped <- make_grouped()
  subgroup_sizes <- rep(4, nrow(grouped))

  r_stats <- qcc::stats.R(grouped, subgroup_sizes)
  r_sd <- qcc::sd.R(grouped, subgroup_sizes)
  r_limits <- qcc::limits.R(
    center = r_stats$center,
    std.dev = r_sd,
    sizes = subgroup_sizes,
    conf = 0.95
  )

  expect_equal(length(r_stats$statistics), nrow(grouped))
  expect_equal(ncol(r_limits), 2)
  expect_true(all(r_limits[, "UCL"] >= r_limits[, "LCL"]))

  max_r_size <- length(qcc::qcc.options("se.R.unscaled")) + 1
  expect_error(
    qcc::limits.R(
      center = r_stats$center,
      std.dev = r_sd,
      sizes = max_r_size,
      nsigmas = 3
    ),
    "group size must be less than"
  )
})
