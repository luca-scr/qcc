theme_qcc <- function(...) {
  bg_margin <- qcc.options("bg.margin")
  bg_figure <- qcc.options("bg.figure")

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
