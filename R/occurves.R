#' Operating Characteristic Function
#'
#' Draws the operating characteristic curves for a `'qcc'` object.
#'
#' An operating characteristic curve graphically provides information about the
#' probability of not detecting a shift in the process. `ocCurves` is a
#' generic function which calls the proper function depending on the type of
#' `'qcc'` object. Further arguments provided through `...` are
#' passed to the specific function depending on the type of chart.
#'
#' The probabilities are based on the conventional assumptions about process
#' distributions: the normal distribution for `"xbar"`, `"R"`, and
#' `"S"`, the binomial distribution for `"p"` and `"np"`, and
#' the Poisson distribution for `"c"` and `"u"`. They are all
#' sensitive to departures from those assumptions, but to varying degrees. The
#' performance of the `"S"` chart, and especially the `"R"` chart,
#' are likely to be seriously affected by longer tails.
#' 
#' @export
#' @param object an object of class `'qcc'`.
#' @param size a vector of values specifying the sample sizes for which to draw
#' the OC curves.
#' @param shift,multiplier a vector of values specifying the shift or
#' multiplier values (in units of sigma).
#' @param nsigmas a numeric value specifying the number of sigmas to use for
#' computing control limits; if `nsigmas` is `NULL`,
#' `object$conf` is used to set up probability limits.
#' @param x an object of class `'ocCurves'`.
#' @param digits the number of significant digits to use.
#' @param what a string specifying the quantity to plot on the y-axis. Possible
#' values are `"beta"` for the probability of not detecting a shift, and
#' `"ARL"` for the average run length.
#' @param title a character string specifying the main title. Set `title =
#' NULL` to remove the title.
#' @param xlab,ylab a string giving the label for the x-axis and the y-axis.
#' @param lty,lwd,col values or vector of values controlling the line type,
#' line width and colour of curves.
#' @param ... catches further ignored arguments.
#' @return The function returns an object of class `'ocCurves'` which
#' contains a matrix or a vector of beta values (the probability of type II
#' error) and ARL (average run length).
#' @author Luca Scrucca
#' @seealso [qcc()]
#' @references `r refs("mason_young_2002", "montgomery2013", "ryan_2011", "scrucca_2004", "wetherill_brown_1991")`
#' @examples
#'
#' data(pistonrings)
#' diameter  = qccGroups(diameter, sample, data = pistonrings)
#' oc  = ocCurves.xbar(qcc(diameter, type="xbar", nsigmas=3))
#' oc
#' plot(oc)
#'
#' data(orangejuice)
#' oc  = with(orangejuice,
#'            ocCurves(qcc(D[trial], sizes=size[trial], type="p")))
#' oc
#' plot(oc)
#'
#' data(circuit)
#' oc  = with(circuit,
#'            ocCurves(qcc(x[trial], sizes=size[trial], type="c")))
#' oc
#' plot(oc)
#'
ocCurves <- function(object, ...)
{
# Compute and draws the operating characteristic curves for a qcc object 

  if ((missing(object)) | (!inherits(object, "qcc")))
     stop("an object of class 'qcc' is required")

  size <- unique(object$sizes)
  if (length(size)>1)
     stop("Operating characteristic curves available only for equal sample sizes!")

  out <- switch(object$type,
                xbar = ocCurves.xbar(object, ...),
                R    = ocCurves.R(object, ...),
                S    = ocCurves.S(object, ...),
                np   =,
                p    = ocCurves.p(object, ...),
                u    =,
                c    = ocCurves.c(object, ...))
  if(is.null(out))
    stop("Operating characteristic curves not available for this type of chart.")
  return(out)
}


# Assemble OC results while preserving chart-specific fields and matrix labels.
.new_oc_curves <- function(type, beta, grid, grid.name, row.label, size = NULL)
{
  colnames(beta) <- if (is.null(size)) "beta" else size
  rownames(beta) <- sprintf(paste0("%.", max(nchar(sub(".*\\.", "", grid))), "f"), grid)
  names(dimnames(beta)) <- c(row.label, if (is.null(size)) "" else "sample size")

  ARL <- 1 / (1 - beta)
  if (is.null(size))
    colnames(ARL) <- "ARL"

  out <- list(type = type)
  if (!is.null(size))
    out$size <- size
  out[[grid.name]] <- grid
  out$beta <- beta
  out$ARL <- ARL
  class(out) <- "ocCurves"
  return(out)
}


#' @rdname ocCurves
#' @export
ocCurves.xbar <- function(object, 
                          size = c(1,5,10,15,20), 
                          shift = seq(0, 5, by = 0.1), 
                          nsigmas = object$nsigmas, ...)
{
# Compute beta and ARL for xbar-chart with nsigmas limits. 

  if (!(object$type == "xbar"))
     stop("not a 'qcc' object of type \"xbar\".")

  n <- unique(object$sizes)
  if (length(n) > 1)
     stop("Operating characteristic curves available only for equal sample sizes!")
  size <- sort(unique(c(n, size)))
  if (is.null(nsigmas))
     nsigmas <- qnorm(1 - (1 - object$confidence.level) / 2)

  beta <- matrix(as.double(NA), nrow = length(shift), ncol = length(size))
  for(i in 1:length(size))
  {
    beta[,i] <- pnorm(nsigmas-shift*sqrt(size[i])) - 
                pnorm(-nsigmas-shift*sqrt(size[i]))
  }
  .new_oc_curves(object$type, beta, shift, "shift", "shift (StdDev)", size = size)
}

#' @rdname ocCurves
#' @export
ocCurves.R <- function(object, 
                       size = c(2,5,10,15,20), 
                       multiplier = seq(1, 6, by = 0.1),
                       nsigmas = object$nsigmas, ...)
{
# Computes the operating-characteristic curves for the R-chart with nsigmas
# limits. The values on the vertical axis give the probability of not detecting
# a change from sigma to c*sigma on the first sample following the change.

  if (!(object$type=="R"))
     stop("not a `qcc' object of type \"R\".")

  n <- unique(object$sizes)
  if (length(n) > 1)
     stop("Operating characteristic curves available only for equal sample sizes!")
  size <- sort(unique(c(n, size)))

  if(is.null(nsigmas))
  { 
    tail.prob <- (1 - object$confidence.level) / 2
    beta.fun1 <- function(c, n, p)
    {
      lcl <- qtukey(p, n, Inf)
      ucl <- qtukey(p, n, Inf, lower.tail = FALSE)
      ptukey(ucl / c, n, Inf) - ptukey(lcl / c, n, Inf)
    }
    beta <- outer(multiplier, size, beta.fun1, tail.prob)
  }
  else
  { 
    beta.fun2 <- function(c, n, conf) {
      lcl <- pmax(0, .d2(n) - conf * .d3(n))
      ucl <- .d2(n) + conf * .d3(n)
      ptukey(ucl / c, n, Inf) - ptukey(lcl / c, n, Inf)
    }
    beta <- outer(multiplier, size, beta.fun2, nsigmas)
  }
  .new_oc_curves(object$type, beta, multiplier, "multiplier", "scale multiplier", size = size)
}

#' @rdname ocCurves
#' @export
ocCurves.S <- function(object, 
                       size = c(2,5,10,15,20), 
                       multiplier = seq(1,6,by=0.1),
                       nsigmas = object$nsigmas, ...)
{
# Computes the operating-characteristic curves for the S-chart with nsigmas
# limits. The values on the vertical axis give the probability of not detecting
# a change from sigma to c*sigma on the first sample following the change.

  if (!(object$type=="S"))
     stop("not a `qcc' object of type \"S\".")

  n <- unique(object$sizes)
  if (length(n) > 1)
     stop("Operating characteristic curves available only for equal sample sizes!")
  size <- sort(unique(c(n, size)))

  if(is.null(nsigmas))
  { 
    tail.prob <- (1 - object$confidence.level) / 2
    beta.fun1 <- function(c, n, p)
    {
      lcl <- sqrt(qchisq(p, n - 1) / (n - 1))
      ucl <- sqrt(qchisq(1 - p, n - 1) / (n - 1))
      pchisq((n-1)*(ucl/c)^2, n-1) - pchisq((n-1)*(lcl/c)^2, n-1)
    }
    beta <- outer(multiplier, size, beta.fun1, tail.prob)
  }
  else
  { 
    beta.fun2 <- function(c, n, nsigmas)
    {
      center <- .c4(n)
      tol <- sqrt(1 - .c4(n)^2)
      lcl <- pmax(0, center - nsigmas * tol)
      ucl <- center + nsigmas * tol
      pchisq((n-1)*(ucl/c)^2, n-1) - pchisq((n-1)*(lcl/c)^2, n-1)
    }
    beta <- outer(multiplier, size, beta.fun2, nsigmas)
  }
  .new_oc_curves(object$type, beta, multiplier, "multiplier", "scale multiplier", size = size)
}

#' @rdname ocCurves
#' @export
ocCurves.p <- function(object, ...)
{
  if (!(object$type=="p" | object$type=="np"))
     stop("not a `qcc' object of type \"p\" or \"np\".")

  size <- unique(object$sizes)
  if (length(size) > 1)
     stop("Operating characteristic curves available only for equal sample sizes!")

  if (is.null(object$limits))
     stop("the `qcc' object does not have control limits!")
  limits <- object$limits
  p <- seq(0, 1, length=101)

  if(object$type=="p") 
  { 
    UCL <- min(floor(size*limits[,2]), size)
    LCL <- max(floor(size*limits[,1]), 0) 
  } else
  { 
    UCL <- min(floor(limits[,2]), size)
    LCL <- max(floor(limits[,1]), 0) 
  }
  beta <- matrix(pbinom(UCL, size, p) - pbinom(LCL-1, size, p), ncol = 1)
  warning("Some computed values for the type II error have been rounded due to the discreteness of the binomial distribution. Thus, some ARL values might be meaningless.")

  .new_oc_curves(object$type, beta, p, "p", "fraction nonconforming")
}

#' @rdname ocCurves
#' @export
ocCurves.c <- function(object, ...)
{
  if (!(object$type=="c" | object$type=="u"))
     stop("not a `qcc' object of type \"c\" or \"u\".")

  size <- unique(object$sizes)
  if (length(size) > 1)
     stop("Operating characteristic curves available only for equal sample size!")

  if (is.null(object$limits))
     stop("the `qcc' object does not have control limits!")
  limits <- object$limits
  CL  <- object$center
  std.dev <- object$std.dev
  if(object$type=="c") 
  { 
    max.lambda <- ceiling(CL+10*std.dev)
    UCL <- floor(limits[1,2])
    LCL <- floor(limits[1,1])
  } else
  { 
    max.lambda <- ceiling(CL*size+10*std.dev*sqrt(size))[1]
    UCL <- floor(size*limits[1,2])
    LCL <- floor(size*limits[1,1])
  }
  lambda <- seq(0, max.lambda)
  beta <- matrix(ppois(UCL, lambda) - ppois(LCL-1, lambda), ncol = 1)

  warning("Some computed values for the type II error have been rounded due to the discreteness of the Poisson distribution. Thus, some ARL values might be meaningless.")
  .new_oc_curves(object$type, beta, lambda, "lambda", "average nonconforming")
}

#' @rdname ocCurves
#' @export
#' @export print.ocCurves
print.ocCurves <- function(x, digits =  getOption("digits"), ...)
{
  object <- x   # Argh.  Really want to use 'object' anyway
  cat(cli::rule(left = crayon::bold("Operating Characteristic Curves"), 
                width = min(getOption("width"),50)), "\n\n")
  cat("Chart type: ", object$type, "\n")
  cat("\nProb. type II error (beta):\n")
  .printShortMatrix(zapsmall(object$beta,digits = 4), head = 3, tail = 2)
  cat("\nAverage run length (ARL):\n")
  .printShortMatrix(zapsmall(object$ARL,digits = 2), head = 3, tail = 2)
  
  invisible()
}


# chart-specific plotting defaults.
.oc_curves_plot_data <- function(object, what)
{
  response <- object[[what]]
  axis <- switch(object$type,
    xbar = list(x = object$shift,
                label = "Process shift (StdDev)",
                scale = list(breaks = unique(as.integer(object$shift)))),
    R =,
    S = list(x = object$multiplier,
             label = "Process scale multiplier",
             scale = list(breaks = unique(as.integer(object$multiplier)))),
    p =,
    np = list(x = object$p,
              label = "Fraction nonconforming",
              scale = list(breaks = seq(0, 1, by = 0.2))),
    c =,
    u = list(x = object$lambda,
             label = "Average nonconforming",
             scale = list(n.breaks = 7)))

  grouped <- object$type %in% c("xbar", "R", "S")
  ncurves <- if(grouped) length(object$size) else 1L
  df <- data.frame(x = rep(axis$x, times = ncurves),
                   y = c(response),
                   size = if(grouped)
                     factor(rep(object$size, each = length(axis$x)))
                   else NA)

  list(data = df, grouped = grouped, ncurves = ncurves,
       xlab = axis$label,
       ylab = if(what == "beta") "Prob. type II error" else "ARL",
       xscale = axis$scale,
       yscale = if(what == "beta") list(breaks = seq(0, 1, by = 0.1))
                else list(n.breaks = 9))
}


#' @rdname ocCurves
#' @export
#' @export plot.ocCurves
plot.ocCurves <- function(x, what = c("beta", "ARL"),
                          title, xlab, ylab, lty, lwd, col,
                          ...)
{
# Draw the operating-characteristic curves. 
# The values on the vertical axis give the probability of not detecting
# a shift of c*sigma in the mean on the first sample following the shift.

  object <- x  # Argh.  Really want to use 'object' anyway
  stopifnot(inherits(object, "ocCurves"))
  what <- match.arg(what, choices = eval(formals(plot.ocCurves)$what), 
                    several.ok = FALSE)
  plotting <- .oc_curves_plot_data(object, what)
  if(missing(title))
    title <- paste("OC curves for", object$type, "chart")
  if(missing(xlab))
    xlab <- plotting$xlab
  if(missing(ylab))
    ylab <- plotting$ylab
  if(missing(lty))
    lty <- rep(1, plotting$ncurves)
  if(missing(lwd))
    lwd <- rep(1, plotting$ncurves)
  if(missing(col))
    col <- blues.colors(plotting$ncurves)

  plot <- ggplot(plotting$data, aes(x = .data[["x"]], y = .data[["y"]]))
  if(plotting$grouped)
  {
    plot <- plot +
      geom_line(aes(linetype = .data[["size"]],
                    linewidth = .data[["size"]],
                    colour = .data[["size"]])) +
      scale_linetype_manual(values = lty) +
      scale_linewidth_manual(values = lwd) +
      scale_colour_manual(values = col) +
      labs(linetype = "Sample size:",
           linewidth = "Sample size:",
           colour = "Sample size:")
  } else
  {
    plot <- plot +
      geom_line(linewidth = lwd[1], col = col[1], lty = lty[1])
  }

  plot <- plot +
    labs(title = title, x = xlab, y = ylab) +
    do.call(scale_x_continuous, plotting$xscale) +
    do.call(scale_y_continuous, plotting$yscale)

  plot <- plot + 
    theme_light() + 
    theme(
      plot.background = element_rect(
        fill = getOption("qcc.bg.margin"),
        color = getOption("qcc.bg.margin")
      ),
      panel.background = element_rect(
        fill = getOption("qcc.bg.figure")
      ),
      plot.title = element_text(face = "bold", size = 11),
      legend.position = c(0.9,0.8),
      legend.text.align = 0,
      axis.text.y = element_text(
        angle = 90, 
        margin = margin(l = 5, r = 5),
        hjust = 0.5, vjust = 0.5)
    )
  
  return(plot)
}
