#' Shared Control Chart (Shewhart, EWMA, CUSUM) Documentation
#'
#' @name chart_common
#' @param data a data frame, a matrix or a vector containing observed data for
#'   the variable to chart. Each row of a data frame or a matrix, and each value
#'   of a vector, refers to a sample or ''rationale group''.
#' @param newdata a data frame, matrix or vector, as for the `data`
#'   argument, providing further data to plot but not included in the
#'   computations.
#' @param newsizes a vector as for the `sizes` argument providing further
#'   data sizes to plot but not included in the computations.
#' @param center a value specifying the center of group statistics or target.
#' @keywords internal
NULL

#' Shared Plotting Documentation
#'
#' @name plot_common
#' @param xtime a vector of date-time values as returned by
#'   [Sys.time()] and [Sys.Date()]. If provided it is used
#'   for x-axis so it must be of the same length as the statistic charted.
#' @param add.stats a logical value indicating whether statistics and other
#'   information should be printed at the bottom of the chart.
#' @param chart.all a logical value indicating whether both statistics for
#'   `data` and for `newdata` (if given) should be plotted.
#' @param fill a logical value specifying if the in-control area should be
#'   filled with the color specified in `qcc.options("zones")$fill`.
#' @param title a character string specifying the main title. Set `title =
#'   NULL` to remove the title.
#' @param xlab a string giving the label for the x-axis.
#' @param ylab a string giving the label for the y-axis.
#' @param xlim a numeric vector specifying the limits for the x-axis.
#' @param ylim a numeric vector specifying the limits for the y-axis.
#' @param digits the number of significant digits to use.
#' @keywords internal
NULL 

# TODO: Select parameters in `@inheritParams plot_common` calls when [https://github.com/r-lib/roxygen2/issues/1879] gets fixed.
# As of writing this, `roxygen2` does not allow multiple filtered `@inheritParams` in the same `.Rd`
# `plot.<object>()` shares same `.Rd` as `<object>()`,
# So we can't filter both `@inheritParams spc_common` and `@inheritParams plot_common` until roxygen/issue-1879 gets fixed.



#' Shared Chart Methods Documentation 
#'
#' @name spc_common
#' @param data The observed data values.
#' @param center The sample/group center statistic.
#' @param sizes The sample sizes.
#' @param std.dev The within-group standard deviation.
#' @param nsigmas Number of sigmas used to compute control limits.
#' Ignored when the `conf` argument is provided.
#' @param conf Confidence level used to compute control limits.
#' Must be a numeric value in \eqn{(0,1)}.
#' @author Luca Scrucca
#' @seealso [qcc()]
#' @keywords internal
NULL


#' Package bibliography
#'
#' A named list of bibliographic references used throughout the package
#' documentation.
#'
#' Reference keys are passed to `refs()` to generate formatted reference lists.
#'
#' @format A named list of character vectors. Each vector contains the
#'   bibliographic components of one reference.
#'
#' @keywords internal
#' @noRd
.bibliography <- c(
  montgomery2013 = paste(
    "Montgomery, D.C. (2013).",
    "*Introduction to Statistical Quality Control*, 7th ed.",
    "New York: John Wiley & Sons."
  ),
  wetherill_brown_1991 = paste(
    "Wetherill, G. B., and Brown, D. W. (1991).",
    "*Statistical Process Control*.",
    "New York: Chapman & Hall."
  ),
  kaminsky_1992 = paste(
    "Kaminsky, F. C., et al. (1992).",
    "*Statistical Control Charts Based on a Geometric Distribution*.",
    "*Journal of Quality Technology*, **24**, 63--69."
  ),
  yang_2002 = paste(
    "Yang, Z., et al. (2002).",
    "On the Performance of Geometric Charts with Estimated Control Limits.",
    "*Journal of Quality Technology*, **34**, 448--458."
  ),
  burr_1969 = paste(
    "Burr, I. W. (1969).",
    "Control charts for measurements with varying sample sizes.",
    "*Journal of Quality Technology*, **1**(3), 163--167."
  ),
  ryan_2011 = paste(
    "Ryan, T. P. (2011).",
    "*Statistical Methods for Quality Improvement*, 3rd ed.",
    "New York: John Wiley & Sons, Inc."
  ),
  mason_young_2002 = paste(
    "Mason, R. L., and Young, J. C. (2002).",
    "*Multivariate Statistical Process Control with Industrial Applications*.",
    "SIAM."
  ),
  scrucca_2004 = paste(
    "Scrucca, L. (2004).",
    "qcc: an R package for quality control charting and statistical process control.",
    "*R News*, **4**(1), 11--17."
  ),
  bissell_1990 = paste(
    "Bissell, A. F. (1990).",
    "How reliable is your capability index?",
    "*Applied Statistics*, **39**, 331--340."
  ),
  boyles_1991 = paste(
    "Boyles, R. A. (1991).",
    "The Taguchi capability index.",
    "*Journal of Quality Technology*, **23**, 107--126."
  ),
  chou_owen_borrego_1990 = paste(
    "Chou, Y., Owen, D. B., and Borrego, S. A. (1990).",
    "Lower confidence limits on process capability indices.",
    "*Journal of Quality Technology*, **22**, 223--229."
  ),
  cano_moguerza_redchuk_2012 = paste(
    "Cano, Emilio L., Moguerza, Javier M. and Redchuk, Andres. 2012.",
    "*Six Sigma with {R}. Statistical Engineering for ProcessImprovement*,",
    "Use R!, vol. 36. Springer, New York."
  )
)

#' Format package references
#'
#' Retrieves selected entries from `.bibliography` and formats them as a
#' Markdown bullet list for inclusion in package documentation.
#'
#' @examples
#' refs("montgomery2013")
#' refs("kaminsky_1992", "yang_2002")
#'
#' @keywords internal
#' @noRd
refs <- function(...) {
  keys <- c(...)

  unknown <- setdiff(keys, names(.bibliography))

  if (length(unknown)) {
    stop(
      "Unknown reference key(s): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  paste0(
    "- ",
    unname(.bibliography[keys]),
    collapse = "\n\n"
  )
}
