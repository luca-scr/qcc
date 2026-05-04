#' Supplementary non-randomness diagnostics for qcc charts
#'
#' Detects patterns 5-8 from Nelson (1984) as a post-hoc analysis on an
#' existing \code{qcc} object. The original object and its \code{$violations}
#' field are never modified.
#'
#' @param object An object of class \code{"qcc"}.
#' @param pattern5 Logical. If \code{TRUE}, detect 15 consecutive points in
#'   Zone C (stratification signal).
#' @param pattern6 Logical. If \code{TRUE}, detect 8 consecutive points
#'   outside Zone C (mixture signal).
#' @param pattern7 Logical. If \code{TRUE}, detect 14 consecutive alternating
#'   points (oscillation signal).
#' @param pattern8 Logical. If \code{TRUE}, detect 7 consecutive trending
#'   points (drift signal).
#'
#' @return An object of class \code{"qcc_diagnostics"} containing:
#'   \item{qcc_object}{the original \code{qcc} object, unmodified}
#'   \item{stats}{full statistics vector (phase I + phase II)}
#'   \item{zone_c}{matrix of Zone C limits (±1σ) per observation}
#'   \item{pattern5}{integer vector of violation indices, or \code{NULL} if
#'     not requested}
#'   \item{pattern6}{integer vector of violation indices, or \code{NULL} if
#'     not requested}
#'   \item{pattern7}{integer vector of violation indices, or \code{NULL} if
#'     not requested}
#'   \item{pattern8}{integer vector of violation indices, or \code{NULL} if
#'     not requested}
#'
#' @references
#' Nelson, L.S. (1984). The Shewhart Control Chart: Tests for Special Causes.
#' \emph{Journal of Quality Technology} 16(4), 237-239.
#'
#' Montgomery, D.C. (2009). \emph{Introduction to Statistical Quality Control},
#' 6th ed. Wiley. Section 5.4.
#'
#' @seealso \code{\link{qcc}}, \code{\link{qccRules}}
#'
#' @export
qcc_diagnostics <- function(object,
                              pattern5 = FALSE,
                              pattern6 = FALSE,
                              pattern7 = FALSE,
                              pattern8 = FALSE)
{
  if (!inherits(object, "qcc"))
    stop("'object' must be of class 'qcc'")

  s      <- c(object$statistics, object$newstats)
  zone_c <- .qcc_zone_c(object)

  result <- list(
    qcc_object = object,
    stats      = s,
    zone_c     = zone_c,
    pattern5   = if (pattern5) .detect_p5(s, zone_c) else NULL,
    pattern6   = if (pattern6) .detect_p6(s, zone_c) else NULL,
    pattern7   = if (pattern7) .detect_p7(s)          else NULL,
    pattern8   = if (pattern8) .detect_p8(s)          else NULL,
    call       = match.call()
  )
  class(result) <- "qcc_diagnostics"
  result
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Compute Zone C limits (±1σ) for all observations using the same
# limits.<type>() infrastructure as the rest of the package.
.qcc_zone_c <- function(object)
{
  all_sizes <- c(object$sizes, object$newsizes)
  limits_fn <- paste0("limits.", object$type)
  if (!exists(limits_fn, mode = "function"))
    stop(paste("function", limits_fn, "is not defined"))

  zone_c <- do.call(limits_fn,
                    list(center  = object$center,
                         std.dev = object$std.dev,
                         sizes   = all_sizes,
                         nsigmas = object$nsigmas * (1 / 3)))

  # When limits are fixed (1 row), expand to full series length
  n <- length(c(object$statistics, object$newstats))
  if (nrow(zone_c) == 1L)
    zone_c <- zone_c[rep(1L, n), , drop = FALSE]

  zone_c
}

# embed()-based window detector: returns indices of the LAST point in every
# window where all run_len values satisfy the condition.
.run_all_true <- function(condition, run_len)
{
  if (length(condition) < run_len) return(integer(0L))
  windows <- embed(condition, run_len)
  as.integer(which(rowSums(windows) == run_len) + (run_len - 1L))
}

# ---------------------------------------------------------------------------
# Pattern detectors
# ---------------------------------------------------------------------------

# Pattern 5: 15 consecutive points in Zone C (stratification)
.detect_p5 <- function(s, zone_c, run_len = 15L)
{
  in_c <- s >= zone_c[, 1] & s <= zone_c[, 2]
  .run_all_true(in_c, run_len)
}

# Pattern 6: 8 consecutive points outside Zone C (mixture)
.detect_p6 <- function(s, zone_c, run_len = 8L)
{
  out_c <- s < zone_c[, 1] | s > zone_c[, 2]
  .run_all_true(out_c, run_len)
}

# Pattern 7: 14 consecutive alternating points (oscillation)
# 14 alternating points require 13 consecutive differences and 12 consecutive
# sign-change pairs in those differences.
.detect_p7 <- function(s, run_len = 14L)
{
  n <- length(s)
  if (n < run_len) return(integer(0L))
  diffs    <- diff(s)
  # TRUE where consecutive differences have opposite signs (product < 0)
  sign_alt <- diffs[-length(diffs)] * diffs[-1] < 0
  n_pairs  <- run_len - 2L   # 12 consecutive sign-change pairs needed
  if (length(sign_alt) < n_pairs) return(integer(0L))
  # .run_all_true returns the last index in sign_alt for each window.
  # sign_alt[j] involves s[j] through s[j+2], so the last point of the
  # 14-point window in s is 2 positions beyond the last sign_alt index.
  as.integer(.run_all_true(sign_alt, n_pairs) + 2L)
}

# Pattern 8: 7 consecutive trending points (drift)
# 7 points in a strictly monotone sequence require 6 same-sign differences.
.detect_p8 <- function(s, run_len = 7L)
{
  n <- length(s)
  if (n < run_len) return(integer(0L))
  diffs <- diff(s)
  steps <- run_len - 1L   # 6 same-sign differences needed
  hits  <- c(.run_all_true(diffs > 0, steps),
             .run_all_true(diffs < 0, steps))
  as.integer(sort(unique(hits)))
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' @export
print.qcc_diagnostics <- function(x, ...)
{
  object <- x
  cat(cli::rule(left = crayon::bold("QCC Supplementary Diagnostics"),
                width = min(getOption("width"), 50L)), "\n\n")

  desc <- c(
    "5" = "15 consecutive in Zone C     (stratification)",
    "6" = "8 consecutive outside Zone C (mixture)",
    "7" = "14 consecutive alternating   (oscillation)",
    "8" = "7 consecutive trending       (drift)"
  )

  for (pid in c("5", "6", "7", "8")) {
    idx   <- object[[paste0("pattern", pid)]]
    label <- paste0("Pattern ", pid, " — ", desc[pid])
    if (is.null(idx)) {
      cat(label, "\n  not requested\n\n")
    } else if (length(idx) == 0L) {
      cat(label, "\n  no violations detected\n\n")
    } else {
      cat(label, "\n")
      cat("  ", length(idx), "violation(s) at index(es):",
          paste(idx, collapse = ", "), "\n\n")
    }
  }
  invisible(object)
}

#' @export
plot.qcc_diagnostics <- function(x,
                                  xtime     = NULL,
                                  add.stats = qcc.options("add.stats"),
                                  title, xlab, ylab, xlim, ylim,
                                  digits    = getOption("digits"), ...)
{
  object  <- x
  qcc_obj <- object$qcc_object
  s       <- object$stats

  # x-axis positions must mirror plot.qcc() logic
  groups <- if (is.null(xtime)) seq_along(s) else xtime

  # Build base chart via plot.qcc
  base_args <- list(x         = qcc_obj,
                    xtime     = xtime,
                    add.stats = add.stats,
                    digits    = digits)
  if (!missing(title)) base_args$title <- title
  if (!missing(xlab))  base_args$xlab  <- xlab
  if (!missing(ylab))  base_args$ylab  <- ylab
  if (!missing(xlim))  base_args$xlim  <- xlim
  if (!missing(ylim))  base_args$ylim  <- ylim
  p <- do.call(plot, base_args)

  # Visual encoding: color + shape per pattern (distinct from WER markers)
  patt_col   <- c("5" = "#E69F00", "6" = "#9370DB",
                   "7" = "#228B22", "8" = "#C0392B")
  patt_shape <- c("5" = 24L, "6" = 23L, "7" = 4L, "8" = 15L)
  patt_label <- c("5" = "P5: 15 in Zone C",
                   "6" = "P6: 8 outside Zone C",
                   "7" = "P7: 14 alternating",
                   "8" = "P8: 7 trending")

  caption_parts <- character(0L)

  for (pid in c("5", "6", "7", "8")) {
    idx <- object[[paste0("pattern", pid)]]
    if (is.null(idx)) next

    n_viol <- length(idx)
    caption_parts <- c(caption_parts,
                        paste0(patt_label[pid], " [", n_viol, "]"))
    if (n_viol == 0L) next

    vdf <- data.frame(x = groups[idx], y = s[idx])
    p[[1]] <- p[[1]] +
      ggplot2::geom_point(
        data    = vdf,
        mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
        colour  = unname(patt_col[pid]),
        shape   = unname(patt_shape[pid]),
        size    = 3,
        stroke  = 1.2
      )
  }

  if (length(caption_parts) > 0L)
    p[[1]] <- p[[1]] +
      ggplot2::labs(caption = paste(caption_parts, collapse = "   "))

  p
}
