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
#' @param ... Additional ignored arguments.
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
