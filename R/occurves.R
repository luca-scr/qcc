#-------------------------------------------------------------------#
#                                                                   #
#          Operating Characteristic Curves                          #
#                                                                   #
#-------------------------------------------------------------------#

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
  colnames(beta) <- size
  rownames(beta) <- sprintf(paste0("%.", max(nchar(sub(".*\\.", "", shift))), "f"), shift)
  names(dimnames(beta)) <- c("shift (StdDev)", "sample size")

  ARL <- 1/(1-beta)

  out <- list(type = object$type, 
              size = size, shift = shift, 
              beta = beta, ARL = ARL) 
  class(out) <- "ocCurves"
  return(out)
}

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
    exp.R.unscaled <- qcc.options("exp.R.unscaled")
    se.R.unscaled <- qcc.options("se.R.unscaled")
    Rtab <- min(length(exp.R.unscaled), length(se.R.unscaled))
    if (any(size > Rtab))
      stop(paste("group size must be less than",
                 Rtab + 1, "when giving nsigmas"))
    beta.fun2 <- function(c, n, conf)
    {
      d2 <- exp.R.unscaled[n]
      d3 <- se.R.unscaled[n]
      lcl <- pmax(0, d2 - conf * d3)
      ucl <- d2 + conf * d3
      ptukey(ucl / c, n, Inf) - ptukey(lcl / c, n, Inf)
    }
    beta <- outer(multiplier, size, beta.fun2, nsigmas)
  }
  colnames(beta) <- size
  rownames(beta) <- sprintf(paste0("%.", max(nchar(sub(".*\\.", "", multiplier))), "f"), multiplier)
  names(dimnames(beta)) <- c("scale multiplier", "sample size")

  ARL <- 1/(1-beta)

  out <- list(type = object$type, 
              size = size, multiplier = multiplier, 
              beta = beta, ARL = ARL) 
  class(out) <- "ocCurves"
  return(out)
}

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
      center <- qcc.c4(n)
      tol <- sqrt(1 - qcc.c4(n)^2)
      lcl <- pmax(0, center - nsigmas * tol)
      ucl <- center + nsigmas * tol
      pchisq((n-1)*(ucl/c)^2, n-1) - pchisq((n-1)*(lcl/c)^2, n-1)
    }
    beta <- outer(multiplier, size, beta.fun2, nsigmas)
  }
  colnames(beta) <- size
  rownames(beta) <- sprintf(paste0("%.", max(nchar(sub(".*\\.", "", multiplier))), "f"), multiplier)
  names(dimnames(beta)) <- c("scale multiplier", "sample size")

  ARL <- 1/(1-beta)

  out <- list(type = object$type, 
              size = size, multiplier = multiplier, 
              beta = beta, ARL = ARL) 
  class(out) <- "ocCurves"
  return(out)
}

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
  rownames(beta) <- sprintf(paste0("%.", max(nchar(sub(".*\\.", "", p))), "f"), p)
  names(dimnames(beta)) <-  c("fraction nonconforming", "\b")
  ARL <- 1/(1-beta)
  
  warning("Some computed values for the type II error have been rounded due to the discreteness of the binomial distribution. Thus, some ARL values might be meaningless.")
  
  out <- list(type = object$type, p = p,
              beta = beta, ARL = ARL) 
  class(out) <- "ocCurves"
  return(out)
}

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
  beta <- ppois(UCL, lambda) - ppois(LCL-1, lambda)
  names(beta) <- sprintf(paste0("%.", max(nchar(sub(".*\\.", "", lambda))), "f"), lambda)
  ARL <- 1/(1-beta)

  warning("Some computed values for the type II error have been rounded due to the discreteness of the Poisson distribution. Thus, some ARL values might be meaningless.")

  out <- list(type = object$type, lambda = lambda,
              beta = beta, ARL = ARL) 
  class(out) <- "ocCurves"
  return(out)
}

print.ocCurves <- function(x, digits =  getOption("digits"), ...)
{
  object <- x   # Argh.  Really want to use 'object' anyway
  cat(cli::rule(left = crayon::bold("Operating Characteristic Curves"), 
                width = min(getOption("width"),50)), "\n\n")
  cat("Chart type                 =", object$type, "\n")
  cat("Prob. type II error (beta) =\n")
  .printShortMatrix(zapsmall(object$beta,digits = 4), head = 3, tail = 2)
  cat("Average run length (ARL)   =\n")
  .printShortMatrix(zapsmall(object$ARL,digits = 2), head = 3, tail = 2)
  
  invisible()
}


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
  if(missing(title))
    title <- paste("OC curves for", object$type, "chart")
  if(missing(ylab))
    ylab <- if(what == "beta") "Prob. type II error" else "ARL"

  if(object$type == "xbar")
  {
    if(missing(xlab))
      xlab <- "Process shift (StdDev)"
    if(missing(lty))
      lty <- rep(1,length(object$size))
    if(missing(lwd))
      lwd <- rep(1,length(object$size))
    if(missing(col))
      col <- blues.colors(length(object$size))
    df <- data.frame(y = c(if(what == "beta") object$beta else object$ARL),
                     size = factor(rep(object$size, 
                                       each = length(object$shift))),
                     shift = rep(object$shift, 
                                 times = length(object$size)))
    plot <- ggplot(df, aes_string(x = "shift", y = "y",
                                  linetype = "size", 
                                  size = "size",
                                  colour = "size")) +
      geom_line() +
      scale_linetype_manual(values = lty) +
      scale_size_manual(values = lwd) +
      scale_colour_manual(values = col) +
      labs(title = title, 
           x = xlab, y = ylab,
           linetype = "Sample size:",
           size = "Sample size:", 
           colour = "Sample size:") +
      scale_x_continuous(breaks = unique(as.integer(object$shift)))
  } else
  if(object$type == "R")
  {
    if(missing(xlab))
      xlab <- "Process scale multiplier"
    if(missing(col))
      col <- blues.colors(length(object$size))
    if(missing(lty))
      lty <- rep(1,length(object$size))
    df <- data.frame(y = c(if(what == "beta") object$beta else object$ARL),
                     size = factor(rep(object$size, 
                                       each = length(object$multiplier))),
                     multiplier = rep(object$multiplier, 
                                      times = length(object$size)))
    plot <- ggplot(df, aes_string(x = "multiplier", 
                                  y = "y",
                                  linetype = "size", 
                                  size = "size",
                                  colour = "size")) +
      geom_line() +
      scale_linetype_manual(values = lty) +
      scale_size_manual(values = lwd) +
      scale_colour_manual(values = col) +
      labs(title = title, 
           x = xlab, y = ylab,
           linetype = "Sample size:",
           size = "Sample size:", 
           colour = "Sample size:") +
      scale_x_continuous(breaks = unique(as.integer(object$multiplier)))
  } else
  if(object$type == "S")
  {
    if(missing(xlab))
      xlab <- "Process scale multiplier"
    if(missing(col))
      col <- blues.colors(length(object$size))
    if(missing(lty))
      lty <- rep(1,length(object$size))
    df <- data.frame(y = c(if(what == "beta") object$beta else object$ARL),
                     size = factor(rep(object$size, 
                                       each = length(object$multiplier))),
                     multiplier = rep(object$multiplier, 
                                      times = length(object$size)))
    plot <- ggplot(df, aes_string(x = "multiplier", 
                                  y = "y", 
                                  linetype = "size", 
                                  size = "size",
                                  colour = "size")) +
      geom_line() +
      scale_linetype_manual(values = lty) +
      scale_size_manual(values = lwd) +
      scale_colour_manual(values = col) +
      labs(title = title, 
           x = xlab, y = ylab,
           linetype = "Sample size:",
           size = "Sample size:", 
           colour = "Sample size:") +
      scale_x_continuous(breaks = unique(as.integer(object$multiplier)))
  } else
  if(object$type == "p" | object$type == "np")
  {
    if(missing(xlab)) xlab <- "Fraction nonconforming"
    if(missing(lty))  lty <- 1
    if(missing(lwd))  lwd <- 1
    if(missing(col))  col <- blues.colors(1)
    df <- data.frame(y = c(if(what == "beta") object$beta else object$ARL),
                     p = object$p)
    plot <- ggplot(df, aes_string(x = "p", y = "y")) +
      geom_line(size = lwd[1], col = col[1], lty = lty[1]) +
      labs(title = title, x = xlab, y = ylab) +
      scale_x_continuous(breaks = seq(0,1,by=0.2)) 
  } else
  if(object$type == "c" | object$type == "u")
  {
    if(missing(xlab)) xlab <- "Average nonconforming"
    if(missing(lty))  lty <- 1
    if(missing(lwd))  lwd <- 1
    if(missing(col))  col <- blues.colors(1)
    df <- data.frame(y = c(if(what == "beta") object$beta else object$ARL),
                     p = object$lambda)
    plot <- ggplot(df, aes_string(x = "p", y = "y")) +
      geom_line(size = lwd[1], col = col[1], lty = lty[1]) +
      labs(title = title, x = xlab, y = ylab) +
      scale_x_continuous(n.breaks = 7)
  }
  
  plot <- plot + 
    if(what == "beta")
      scale_y_continuous(breaks = seq(0,1,by=0.1))
    else
      scale_y_continuous(n.breaks = 9)

  plot <- plot + 
    theme_light() + 
    theme(plot.background = element_rect(fill = qcc.options("bg.margin"),
                                         color = qcc.options("bg.margin")),
          panel.background = element_rect(fill = qcc.options("bg.figure")),
          plot.title = element_text(face = "bold", size = 11),
          legend.position = c(0.9,0.8),
          legend.text.align = 0)
  
  return(plot)
}

