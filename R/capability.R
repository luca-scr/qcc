#' Process capability analysis
#'
#' Computes process capability indices for a `'qcc'` object of type
#' `"xbar"` and plot the histogram.
#'
#' This function calculates confidence limits for \eqn{C_p}{C_p} using the
#' method described by Chou et al. (1990). Approximate confidence limits for
#' \eqn{C_{pl}}{C_pl}, \eqn{C_{pu}}{C_pu}, and \eqn{C_{pk}}{C_pk} are computed
#' using the method in Bissell (1990). Confidence limits for
#' \eqn{C_{pm}}{C_pm} are based on the method of Boyles (1991); this method is
#' approximate and it assumes that the target is midway between the
#' specification limits.
#' 
#' @export
#' @param object a `'qcc'` object of type `"xbar"`
#' @param spec.limits a two-values vector specifying the lower and upper
#' specification limits. For one-sided specification limits, the value of the
#' missing limit must be set to `NA`.
#' @param target a value specifying the target of the process. If missing the
#' value from the `'qcc'` object is used if not `NULL`, otherwise the
#' target is set at the middle value between specification limits.
#' @param std.dev a value specifying the within-group standard deviation. If
#' not provided is taken from the `'qcc'` object.
#' @param nsigmas a numeric value specifying the number of sigmas to use. If
#' not provided is taken from the `'qcc'` object.
#' @param confidence.level a numeric value between 0 and 1 specifying the level
#' to use for computing confidence intervals.
#' @param x an object of class `'processCapability'`.
#' @param add.stats a logical value indicating whether statistics and
#' capability indices should be added at the bottom of the chart.
#' @param breaks a value or a function used to select the number of bins in a
#' histogram. See the help for [nclass.scott()] for more details.
#' @param fill,color values specifying the colour of the filled area and the
#' border used for drawing the histogram.
#' @param title a character string specifying the plot title. Set `title =
#' NULL` to remove the title.
#' @param xlab a character string specifying the label for the x-axis.
#' @param digits the number of significant digits to use.
#' @param ... catches further ignored arguments.
#' @return Invisibly returns a list with components:
#' - `nobs`: number of non-missing observations.
#' - `center`: center.
#' - `std.dev`: standard deviation.
#' - `target`: target.
#' - `spec.limits`: a vector of values giving the lower specification limit
#'   (LSL) and the upper specification limit (USL).
#' - `indices`: a matrix of capability indices (\eqn{C_p}{C_p},
#'   \eqn{C_{pl}}{C_pl}, \eqn{C_{pu}}{C_pu}, \eqn{C_{pk}}{C_pk},
#'   \eqn{C_{pm}}{C_pm}) and the corresponding confidence limits.
#' - `exp`: a vector of values giving the expected fraction, based on a normal
#'   approximation, of the observations less than LSL and greater than USL.
#' - `obs`: a vector of values giving the fraction of observations less than
#'   LSL and greater than USL.
#' @author Luca Scrucca
#' @seealso [qcc()]
#' @references `r refs("bissell_1990", "boyles_1991", "chou_owen_borrego_1990", "montgomery2013", "wetherill_brown_1991")`
#' @examples
#'
#' data(pistonrings)
#' diameter  = qccGroups(data = pistonrings, diameter, sample)
#' q  = qcc(diameter[1:25,], type="xbar", nsigmas=3)
#' pc  = processCapability(q, spec.limits=c(73.95,74.05))
#' pc
#' plot(pc)
#' plot(processCapability(q, spec.limits=c(73.95,74.05), target=74.02))
#' plot(processCapability(q, spec.limits=c(73.99,74.01)))
#' plot(processCapability(q, spec.limits = c(73.99, 74.1)))
#'
processCapability <- function(object, spec.limits, target, 
                              std.dev, nsigmas, 
                              confidence.level = 0.95, ...)
{
# TODO: implement Cpkm

  if ((missing(object)) | (!inherits(object, "qcc")))
     stop("an object of class 'qcc' is required")
  if (!(object$type=="xbar" | object$type=="xbar.one"))
     stop("Process Capability Analysis only available for charts type \"xbar\" and \"xbar.one\" charts")

  x <- as.vector(object$data)
  x <- x[!is.na(x)]
  center <- object$center
  if(missing(std.dev))
    std.dev <- object$std.dev
  overall.std.dev <- stats::sd(x)
  n <- length(x)

  if (missing(spec.limits))
     stop("specification limits must be provided")
  spec.limits <- as.vector(spec.limits)[1:2]
  LSL <- spec.limits[1]
  if(!(is.numeric(LSL) & is.finite(LSL))) LSL <- NA
  USL <- spec.limits[2]
  if(!(is.numeric(USL) & is.finite(USL))) USL <- NA
  if(is.na(LSL) & is.na(USL))
     stop("invalid specification limits")

  # TODO: refactor
  has.target <- (!missing(target))
  if(!has.target) {
    target <- mean(spec.limits, na.rm=TRUE) # FIX: removing NA means target == spec limits
    if(!is.na(LSL) & !is.na(USL)) has.target <- TRUE
    message("target value not provided; using midpoint of specification limits. Cpm and Ppm may be optimistic")
    # TODO: Explain this message in more detail in ?processCapability
  }
     
  if (is.na(LSL))
     { if (target > USL)
           warning("target value larger than one-sided specification limit...") }
  if (is.na(USL))
     { if (target < LSL)
           warning("target value smaller than one-sided specification limit...") }
  if (!is.na(LSL) & !is.na(USL))
     { if (target < LSL || target > USL)
       warning("target value is not within specification limits...") }
       
  if (missing(nsigmas))
     if (is.null(object$nsigmas))
        stop("nsigmas not available in the 'qcc' object. Please provide nsigmas.") 
     else  nsigmas <- object$nsigmas
  
  if (confidence.level < 0 | confidence.level > 1)
     stop("the argument confidence.level must be a value between 0 and 1") 

  # computes process capability indices
  Cp <- (USL - LSL) / (2*nsigmas*std.dev)
  Cp.u <- (USL-center)/(nsigmas*std.dev)
  Cp.l <- (center-LSL)/(nsigmas*std.dev)
  Cp.k <- min(Cp.u, Cp.l, na.rm = TRUE)
  # Cpm <- (USL - LSL) / (2*nsigmas*sqrt(sum((x-target)^2)/(n-1)))
  Cpm <- Cp / sqrt(1+((center-target)/std.dev)^2)
  Pp <- (USL - LSL) / (2*nsigmas*overall.std.dev)
  Pp.u <- (USL-center)/(nsigmas*overall.std.dev)
  Pp.l <- (center-LSL)/(nsigmas*overall.std.dev)
  Pp.k <- min(Pp.u, Pp.l, na.rm = TRUE)
  Ppm <- Pp / sqrt(1+((center-target)/overall.std.dev)^2)

  # compute confidence limits 
  alpha <- 1-confidence.level
  Cp.limits   <- .chisq_limits_cp_family(Cp, n - 1, alpha)
  Cp.u.limits <- .wald_limits_cpk_family(Cp.u, qnorm(confidence.level), n)
  Cp.l.limits <- .wald_limits_cpk_family(Cp.l, qnorm(confidence.level), n)
  Cp.k.limits <- .wald_limits_cpk_family(Cp.k, qnorm(1 - alpha / 2), n)
  df <- n * (1 + ((center - target) / std.dev)^2) / (1 + 2 * ((center - target) / std.dev)^2)
  Cpm.limits <- .chisq_limits_cp_family(Cpm, df, alpha)
  Pp.limits <- .chisq_limits_cp_family(Pp, n - 1, alpha)
  Pp.u.limits <- .wald_limits_cpk_family(Pp.u, qnorm(confidence.level), n)
  Pp.l.limits <- .wald_limits_cpk_family(Pp.l, qnorm(confidence.level), n)
  Pp.k.limits <- .wald_limits_cpk_family(Pp.k, qnorm(1 - alpha / 2), n)
  overall.df <- n * (1 + ((center - target) / overall.std.dev)^2) /
    (1 + 2 * ((center - target) / overall.std.dev)^2)
  Ppm.limits <- .chisq_limits_cp_family(Ppm, overall.df, alpha)

  # limit.names <- (c(alpha/2, 1-alpha/2) * 100) |> round(1) |> paste0("%") # Remove `round`?
  limit.names <- c(paste(round(100*alpha/2, 1), "%", sep=""),
                   paste(round(100*(1-alpha/2), 1), "%", sep=""))
  names(Cp.limits) <- names(Cp.u.limits) <- names(Cp.l.limits) <- names(Cp.k.limits) <-
    names(Cpm.limits) <- names(Pp.limits) <- names(Pp.u.limits) <- names(Pp.l.limits) <-
    names(Pp.k.limits) <- names(Ppm.limits) <- limit.names

  if(is.na(LSL))  exp.LSL <- NA
  else { exp.LSL <- pnorm((LSL-center)/std.dev) * 100
         if(exp.LSL < 0.01) exp.LSL <- 0 }
  if(is.na(USL))  exp.USL <- NA
  else { exp.USL <- (1-pnorm((USL-center)/std.dev)) * 100
         if(exp.USL < 0.01) exp.USL <- 0 }
  obs.LSL <- sum(x<LSL)/n * 100
  obs.USL <- sum(x>USL)/n * 100
  
  tab <- cbind(c(Cp, Cp.l, Cp.u, Cp.k, Cpm, Pp, Pp.l, Pp.u, Pp.k, Ppm),
               rbind(Cp.limits, Cp.l.limits, Cp.u.limits, 
                     Cp.k.limits, Cpm.limits, Pp.limits, Pp.l.limits, Pp.u.limits,
                     Pp.k.limits, Ppm.limits))
  rownames(tab) <- c("Cp", "Cp_l", "Cp_u", "Cp_k", "Cpm",
                     "Pp", "Pp_l", "Pp_u", "Pp_k", "Ppm")
  colnames(tab) <- c("Value", names(Cp.limits))

  out <- list(data = x, data.name = object$data.name,
              center = center, std.dev = std.dev, overall.std.dev = overall.std.dev,
              has.target = has.target, target = target, 
              spec.limits = { sl <- c(LSL, USL)
                              names(sl) <- c("LSL", "USL")
                              sl },
              indices = tab, 
              exp = { exp <- c(exp.LSL, exp.USL)
                      names(exp) <- c("Exp < LSL", "Exp > USL")
                      exp }, 
              obs = { obs <- c(obs.LSL, obs.USL)
                      names(obs) <- c("Obs < LSL", "Obs > USL")
                      obs },
              nobs = n )
  class(out) <- "processCapability"
  return(out)
}

#' @rdname processCapability
#' @method print processCapability
#' @export
#' @export print.processCapability
print.processCapability <- function(x, digits = getOption("digits"), ...)
{
  object <- x   # Argh.  Really want to use 'object' anyway
  cat(cli::rule(left = crayon::bold("Process Capability Analysis"), 
                width = min(getOption("width"),50)), "\n\n")
  
  cat(paste(formatC("Number of obs = ", width=16),
            formatC(length(object$data), width=12, flag="-"),
            formatC("Target = ", width=10),
            ifelse(object$has.target, 
                   formatC(object$target, digits=digits, flag="-"), ""),
            "\n", sep=""))
  cat(paste(formatC("Center        = ", width=16),
            formatC(object$center, digits=digits, width=12, flag="-"),
            formatC("LSL    = ", width=10),
            ifelse(is.na(object$spec.limits[1]), "",
                   formatC(object$spec.limits[1], digits=digits, flag="-")),
            "\n", sep=""))
  cat(paste(formatC("StdDev        = ", width=16),
            formatC(object$std.dev, digits=digits, width=12, flag="-"),
            formatC("USL    = ", width=10),
            ifelse(is.na(object$spec.limits[2]), "",
                   formatC(object$spec.limits[2], digits=digits, flag="-")),
            "\n", sep=""))
  cat(paste(formatC("Overall SD    = ", width=16),
            formatC(object$overall.std.dev, digits=digits, width=12, flag="-"),
            "\n", sep=""))
            
  indices <- object$indices
  names(dimnames(indices)) <- c("Capability indices", "")
  print(indices, digits = 3, na.print="", print.gap=2)
  
  cat("\n")
  cat(paste("Exp<LSL", ifelse(is.na(object$exp[1]), "\t", 
                              paste(signif(object$exp[1], digits=2), "%\t", sep="")), 
            "Obs<LSL", ifelse(is.na(object$obs[1]), "", 
                              paste(signif(object$obs[1], digits=2), "%\n", sep=""))))
  cat(paste("Exp>USL", ifelse(is.na(object$exp[2]), "\t", 
                              paste(signif(object$exp[2], digits=2), "%\t", sep="")),
            "Obs>USL", ifelse(is.na(object$obs[2]), "", 
                              paste(signif(object$obs[2], digits=2), "%\n", sep=""))))
  
  invisible()
}

#' @rdname processCapability
#' @method summary processCapability
#' @export
#' @export summary.processCapability
summary.processCapability <- function(object, ...) 
  print.processCapability(object, ...)

#' @rdname processCapability
#' @method plot processCapability
#' @export
#' @export plot.processCapability
plot.processCapability <- function(x, 
                                   add.stats = getOption("qcc.add.stats"),
                                   breaks = nclass.hist, 
                                   fill = adjustcolor(getOption("qcc.zones")$fill, alpha.f = 0.5), # HACK: too much code for an argument.
                                   color = "white",
                                   title, xlab,
                                   digits = getOption("digits"),
                                   ...)
{
  object <- x  # Argh.  Really want to use 'object' anyway
   if ((missing(object)) | (!inherits(object, "processCapability")))
     stop("an object of class `processCapability' is required")

  nobs <- length(object$data)
  indices <- object$indices[, 1]
  if(is.function(breaks))
    breaks <- breaks(object$data)
  breaks <- as.integer(breaks)
  h <- hist(object$data, breaks = breaks, plot=FALSE)
  xlim <- range(c(h$breaks, object$spec.limits, object$target), na.rm = TRUE)
  xlim <- extendrange(r = xlim, f = 0.1)
  x  <- seq(min(xlim), max(xlim), length=250)
  dx <- dnorm(x, object$center, object$std.dev)
  ylim <- extendrange(c(h$density, dx))
  xlim <- range(c(h$breaks, x))
  
  if(missing(title))
    title <- "Process capability analysis"

  plot <- ggplot() +
    geom_histogram(data = data.frame(data = object$data),
                   aes(x = .data[["data"]], y = after_stat(density)),
                   stat = "bin", breaks = h$breaks,
                   fill = fill, color = color) +
    geom_line(data = data.frame(x, dx), 
              aes(x = x, y = dx)) +
    labs(title = title, subtitle = "", y = "",
         x = if(missing(xlab)) object$data.name else xlab) +
    coord_cartesian(xlim = xlim, ylim = ylim,
                    expand = FALSE, clip = "off") +
    theme_qcc(plot.margin = margin(5, 5, 5, 5)) 
  
  plot <- plot +
    geom_vline(xintercept = object$spec.limits, lty = 2) +
    annotate("text", 
             x = object$spec.limits[1], 
             y = ylim[2],
             label = "LSL", 
             hjust = 0.5, vjust = -0.5, size = 10 * 5/14) +
    annotate("text", 
             x = object$spec.limits[2], 
             y = ylim[2],
             label = "USL", 
             hjust = 0.5, vjust = -0.5, size = 10 * 5/14) 
  
  if(object$has.target)
  { 
    plot <- plot + 
      geom_vline(xintercept = object$target, lty = 3) +
      annotate("text", 
               x = object$target, 
               y = ylim[2],
               label = "Target", 
               hjust = 0.5, vjust = -0.5, size = 10 * 5/14) 
  }
  
  
  if(add.stats) 
  { 
    display <- function(x, digits, suffix = "")
      ifelse(is.na(x), "", paste0(signif(x, digits), suffix))

    sections <- list(
      `Process Summary` = c(
        "Observations" = nobs,
        "Center" = signif(object$center, digits),
        "StdDev" = signif(object$std.dev, digits),
        "Overall SD" = signif(object$overall.std.dev, digits)
      ),
      Specifications = c(
        "Target" = if (object$has.target) signif(object$target, digits) else "",
        "LSL" = signif(object$spec.limits[[1]], digits),
        "USL" = signif(object$spec.limits[[2]], digits)
      ),
      Capability = setNames(
        display(indices[c("Cp", "Cp_l", "Cp_u", "Cp_k", "Cpm")], digits = 3),
        c("C[p]", "C[p*l]", "C[p*u]", "C[p*k]", "C[p*m]") # for geom_text(..., parse = TRUE)
      ),
      Performance = setNames(
        display(indices[c("Pp", "Pp_l", "Pp_u", "Pp_k", "Ppm")], digits = 3),
        c("P[p]", "P[p*l]", "P[p*u]", "P[p*k]", "P[p*m]") # for geom_text(..., parse = TRUE)
      ),
      `Non-conformance` = display(
        c(
          "Exp < LSL" = object$exp[[1]],
          "Exp > USL" = object$exp[[2]],
          "Obs < LSL" = object$obs[[1]],
          "Obs > USL" = object$obs[[2]]
        ),
        digits = 2,
        suffix = "%"
      )
    )

    panels <- chart_footer(
      sections,
      parse = names(sections) %in% c("Capability", "Performance")
    )
    plot <- .add_footer(
      plot,
      panels,
      heights = c(0.83, 0.17),
      widths = c(0.24, 0.16, 0.18, 0.18, 0.24)
    )
  }

  return(plot)
}
  
#' Wald confidence limits for Cpk-family indices
#'
#' Computes approximate two-sided confidence limits for \eqn{C_{pu}}{C_pu},
#' \eqn{C_{pl}}{C_pl}, and \eqn{C_{pk}}{C_pk} using the Wald method described
#' by Bissell (1990).
#'
#' @param idx A numeric scalar giving the point estimate of the capability
#'   index.
#' @param z A numeric scalar giving the normal quantile.
#' @param n A numeric scalar giving the sample size.
#'
#' @return A numeric vector of length two containing the lower and upper
#'   confidence limits. Returns two `NA` values when `idx` is `NA`.
#'
#' @references `r refs("bissell_1990")`
#' @keywords internal
#' @noRd
.wald_limits_cpk_family <- function(idx, z, n)
{
  if (is.na(idx)) return(c(NA_real_, NA_real_))
  idx * (1 + c(-1, 1) * z * sqrt(1 / (9 * n * idx^2) + 1 / (2 * (n - 1))))
}

#' Chi-squared confidence limits for Cp-family indices
#'
#' Computes two-sided confidence limits for \eqn{C_p}{C_p} using the method of
#' Chou et al. (1990), and approximate limits for \eqn{C_{pm}}{C_pm} using the
#' method of Boyles (1991).
#'
#' @param idx A numeric scalar giving the point estimate of the capability
#'   index.
#' @param df A numeric scalar giving the degrees of freedom.
#' @param alpha A numeric scalar giving the total tail probability.
#'
#' @return A numeric vector of length two containing the lower and upper
#'   confidence limits. Returns two `NA` values when `idx` is `NA`.
#'
#' @references `r refs("boyles_1991", "chou_owen_borrego_1990")`
#' @keywords internal
#' @noRd
.chisq_limits_cp_family <- function(idx, df, alpha)
{
  if (is.na(idx)) return(c(NA_real_, NA_real_))
  idx * sqrt(qchisq(c(alpha / 2, 1 - alpha / 2), df) / df)
}
