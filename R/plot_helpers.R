.qcc_default_rules <- function() {
  list(
    col = c(
      "#F03B20", "#EE7600", "#FD8D3C", "#CD5555",
      "#7B3294", "#008837", "#B0BC17", "#C51B7D"
    ),
    pch = c(19, 15, 17, 8, 19, 15, 17, 18)
  )
}

.qcc_default_zones <- function() {
  list(
    fill = "#5E81AC",
    lty = c(2, 2, 2),
    col = grDevices::grey(c(0.1, 0.4, 0.7))
  )
}

.qcc_default_options <- function() {
  list(
    qcc.add.stats = TRUE,
    qcc.chart.all = TRUE,
    qcc.fill = TRUE,
    qcc.rules = .qcc_default_rules(),
    qcc.zones = .qcc_default_zones(),
    qcc.bg.margin = "#EFF0F2",
    qcc.bg.figure = "white",
    qcc.cex = 1,
    qcc.font.stats = 1,
    qcc.cex.stats = 0.9
  )
}

theme_qcc <- function(...) {
  bg_margin <- getOption("qcc.bg.margin")
  bg_figure <- getOption("qcc.bg.figure")

  theme_light() +
    theme(
      plot.background = element_rect(fill = bg_margin, color = bg_margin),
      panel.background = element_rect(fill = bg_figure),
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "none",
      plot.margin = margin(5, 30, 5, 5),
      axis.text.y = element_text(
        angle = 90,
        margin = margin(l = 5, r = 5),
        hjust = 0.5,
        vjust = 0.5
      )
    ) +
    theme(...)
}

theme_qcc_void <- function(...) {
  bg_margin <- getOption("qcc.bg.margin")

  theme_void() +
    theme(
      plot.background = element_rect(fill = bg_margin, color = bg_margin),
      plot.margin = margin(0.5, 0, 0.5, 0, unit = "lines")
    ) +
    theme(...)
}

#' Build a Plot Index for Control Charts
#'
#' Validates plotting coordinates and maps chart rows to displayed groups and
#' phases.
#'
#' @param n_phase1 Number of Phase I rows.
#' @param n_phase2 Number of Phase II rows.
#' @param xtime Optional vector of displayed x-axis coordinates.
#' @param chart_all Logical; whether to retain both phases when Phase II rows
#'   are present.
#'
#' @return A data frame with the positional row, displayed group, and phase.
#'
#' @keywords internal
#' @noRd
qcc_plot_index <- function(n_phase1, n_phase2 = 0L, xtime = NULL, chart_all = TRUE) {
  n <- n_phase1 + n_phase2

  if (!is.null(xtime) &&
      !inherits(xtime, c("numeric", "integer", "Date", "POSIXct", "POSIXt"))) {
    stop(
      "xtime must be of class 'numeric', 'integer', 'Date', 'POSIXct' or 'POSIXt'",
      call. = FALSE
    )
  }

  group <- xtime %||% seq_len(n)
  if (length(group) != n) {
    stop(sprintf("xtime must have length %d", n), call. = FALSE)
  }

  index <- data.frame(
    row = seq_len(n),
    group = group,
    phase = rep(c(1L, 2L), c(n_phase1, n_phase2))
  )

  if (!chart_all && n_phase2 > 0L) {
    index <- index[index$row > n_phase1, , drop = FALSE]
    row.names(index) <- NULL
  }

  index
}

#' Choose an X-Axis Scale for Control Charts
#'
#' Selects a continuous, date, or datetime scale for plotting coordinates.
#'
#' @param x A vector of displayed x-axis coordinates.
#' @param limits A length-two vector giving the x-axis limits.
#' @param n The number of intervals used to compute axis breaks.
#'
#' @return A ggplot2 continuous, date, or datetime position scale.
#'
#' @keywords internal
#' @noRd
scale_x_qcc <- function(x, limits, n = 7L) {
  if (is.numeric(x) || is.integer(x)) {
    scale_x_continuous(breaks = pretty(limits, n = n))
  } else if (inherits(x, "Date")) {
    scale_x_date(breaks = pretty(limits, n = n))
  } else {
    scale_x_datetime(breaks = pretty(limits, n = n))
  }
}

# PERF: Use one ggplot2 object for all footer data, only dowside is changing the user-facing API.
# TODO: figure out a way to format a matrix, in case we add confidence intervals to the footer like JMP.
.add_footer <- function(plot, panels, widths, heights) {
  footer <- patchwork::wrap_plots(plotlist = panels, nrow = 1, widths = widths)
  # TODO: calculate heights from panels$n_rows
  (plot / footer) + patchwork::plot_layout(heights = heights)
}

chart_footer <- function(sections, parse = FALSE) {
  n_rows <- max(lengths(sections))
  row_spacing <- 0.75

  panels <- Map(\(values, section, parse) {
    value_names <- names(values)
    labels <- sprintf("%s:", value_names)
    value_x <- max(nchar(labels, type = "width")) + 1L
    values <- as.character(unname(values))
    data <- data.frame(
      row = n_rows - (seq_along(values) - 1L) * row_spacing,
      label = if (parse) sprintf('%s * ":"', value_names) else labels,
      value = values
    )

    ggplot(data) +
      geom_text(
        aes(y = .data[["row"]], label = .data[["label"]]),
        x = 0, hjust = 0, parse = parse, size = 9, size.unit = "pt"
      ) +
      geom_text(
        aes(y = .data[["row"]], label = .data[["value"]]),
        x = value_x, hjust = 0, size = 9, size.unit = "pt"
      ) +
      labs(title = section) +
      scale_x_continuous(
        limits = c(0, value_x + max(nchar(values, type = "width"))),
        expand = expansion(mult = 0.02)
      ) +
      scale_y_continuous(limits = c(0.5, n_rows + 0.5), expand = expansion(mult = 0)) +
      theme_qcc_void(plot.title = element_text(size = 9, face = "bold"))
  }, sections, names(sections), parse)

  structure(
    panels,
    nrows = n_rows,
    npanels = length(sections),
    class = "footer_panels"
  )
}
