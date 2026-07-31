### Constants of statistical importance in SPC and QC
#
# This file defines d2, d3, and c4.
#
# montgomery 8th Ed also defines A, A3, B3, B4, B5, B6 in Appendix 12.
#  other factors: A2, D1, D2, D3, D4 are tabulated but not defined AFAIK.
#  I think they mostly cause confusion.
#  Strangely, Montgomery does not define the c4 formula.
#
# Minitab support docs additionally defines polynomial approximations for d2, d3, d4,
#  a closed-form expression for c5 and tabulated values for c4' (Not the same as c4).
#
#  I recommend that more constants are added to this file is they can be defined as
#   moments of statistical distributions or have potential for reuse in the package.
#
#  ACCURACY:
#    Numerical integration is more accurate (verified informally) than rational function 
#    approximations suggested by Wardell2025 & Minitab. Thier suggestions probably
#    targets Excel users, not systems with [integrate()] & [ptukey()] equivalents (like R).
#    The only way (I could think of) accuracy could be improved is through tweaking
#    [integrate()] accuracy or short-circuiting to analytic solutions were analytic
#    solutions are defined (as in Wardell2025), That would be the equivalent of
#    a table-lookup.
#
# TODO: Cite Wardell2025 (for d2 and d3) although we don't use their approximation.
# TODO: Cite the original publication for the numerical integration method used in .d2 and .d3
# TODO: Cite the original publication for the closed-form expression used in .c4
# PERFORMANCE: expose a `_d2` argument to `.d3()` and don't use it in `d3()`.
#   this way, `.d3()` can skip recomputing .d2 in formulas that compute both.
#   This would probably be overengineered.
# TODO: d4 -> MMR estimator.
# TODO: c5 -> MVLUE-SD and S-chart SE
# TODO: c4_prime() -> MSSD estimator, due to autocorrelation, c4 would not be 100% correct

#' The \eqn{d_2}{d2} Constant
#'
#' Calculate the \eqn{d_2}{d2} constants for a vector `n` of sample sizes.
#'
#' \eqn{d_2(n)}{d2(n)} is the expected value of the range of n observations
#' from a standard normal distribution.
#' \deqn{d_2(n) = \frac{E(r)}{\sigma}}{d2(n) = E(r)/sigma}
#'
#' @param n A vector of sample sizes.
#' @return A vector of calculated \eqn{d_2}{d2} constants.
#' @references `r refs("cano_moguerza_redchuk_2012")`
#' @family constants for Shewhart charts
#' @seealso For other implementations in R:
#'   [SixSigma::ss.cc.getd2()], [IQCC::d2()], [shewhartr::shewhart_constants()]
#'
#' @export
d2 <- function(n) assert_n(n) |> .d2()

# Analytic solutions for `n` in [2, 5] in Wardell2025
.d2 <- \(n) integrate_ok(
  \(x, n_i) 1 - ptukey(x, n_i, Inf),
  0, Inf, n
)


#' The \eqn{d_3}{d3} Constant
#'
#' Calculate the \eqn{d_3}{d3} constants for a vector `n` of sample sizes.
#'
#' \eqn{d_3(n)}{d3(n)} is the standard deviation of the range of `n` observations
#' from a standard normal distribution.
#' \deqn{d_3(n) = \frac{Stdev(r)}{\sigma}}{d3(n) = stdev(r)/sigma}
#'
#' Not to be confused with the \eqn{D_3}{D3}, which is
#' \deqn{D_4=1-3\frac{d_3}{d_2}}{D4 = 1 - 3 * (d3 / d2)}
#'
#' @param n A vector of sample sizes.
#' @return A vector of calculated \eqn{d_3}{d3} constants.
#' @references `r refs("cano_moguerza_redchuk_2012")`
#' @family constants for Shewhart charts
#' @seealso For other implementations in R:
#'   [SixSigma::ss.cc.getd3()], [IQCC::d3()], [shewhartr::shewhart_constants()]
#'
#' @export
d3 <- function(n) assert_n(n) |> .d3()

# Analytic solutions for `n` in [2, 5] in Wardell2025
.d3 <- \(n) {
  sqrt(
    2 * integrate_ok(
      \(x, n_i) x * (1 - ptukey(x, n_i, Inf)),
      0, Inf, n
    ) - .d2(n)^2
  )
}


#' The \eqn{c_4}{\code{c4}} Constant
#'
#' Calculates the bias-correction factor for the sample standard deviation.
#'
#' For a normally distributed sample of size \eqn{n}, the sample standard
#' deviation \eqn{S} is biased downward,
#' with \deqn{E(S) = c_4(n)\sigma.}{E(s) = c4(n) * sigma}
#' Consequently, \eqn{S / c_4(n)}{S/c4(n)} is an unbiased estimator of the population
#' standard deviation.
#'
#' @param n A vector of sample sizes.
#' @return A vector of calculated \eqn{c_4}{c4} constants.
#' @family constants for Shewhart charts
#' @seealso For other implementations in R:
#'   [SixSigma::ss.cc.getc4()], [IQCC::c4()], [shewhartr::shewhart_constants()]
#' @export
c4 <- function(n) assert_n(n) |> .c4()

# We use  `exp(lgamma(n/2) - lgamma((n - 1)/2))`
# and not `((gamma(n/2))/(gamma((n - 1)/2)))`
# because [gamma()] reteurns `Inf` for n > 171 (On my machine).
.c4 <- \(n) sqrt(2 / (n - 1)) * exp(lgamma(n / 2) - lgamma((n - 1) / 2))


#' Vectorized `integrate()` Wrapper
#'
#' Applies `integrate()` over a parameter vector, preserves missing values, skips
#' repeated parameters, and warns when the estimated absolute error exceeds `max_error`.
#'
#' @keywords internal
#' @noRd
integrate_ok <- \(f, lower, upper, parameter, ..., max_error = 1e-3) {
  parameter_unique <- unique(parameter)

  values <- vapply(
    parameter_unique,
    \(parameter_i) {
      if (is.na(parameter_i))
        return(NA_real_)

      out <- integrate(f, lower, upper, parameter_i, ...)

      if (out$abs.error > max_error) {
        cli_warn(c(
          "Estimated integration error exceeds the threshold.",
          "i" = "Parameter: {.val {parameter_i}}",
          "i" = "Estimated error: {.val {out$abs.error}}",
          "i" = "Threshold: {.val {max_error}}"
        ))
      }

      out$value
    },
    numeric(1)
  )

  values[match(parameter, parameter_unique)]
}
