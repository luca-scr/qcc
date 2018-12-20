#-------------------------------------------------------------------#
#                                                                   #
#          Operating Characteristic Curves                          #
#                                                                   #
#-------------------------------------------------------------------#

ocCurves <- function(object, ...)
{
# Draws the operating characteristic curves for the qcc object 

  if ((missing(object)) | (!inherits(object, "qcc")))
     stop("an object of class 'qcc' is required")

  size <- unique(object$sizes)
  if (length(size)>1)
     stop("Operating characteristic curves available only for equal sample sizes!")

  beta <- switch(object$type,
                 xbar = ocCurves.xbar(object, ...),
                 R    = ocCurves.R(object, ...),
                 S    = ocCurves.S(object, ...),
                 np   =,
                 p    = ocCurves.p(object, ...),
                 u    =,
                 c    = ocCurves.c(object, ...))
  if (is.null(beta))
     stop("Operating characteristic curves not available for this type of chart.")

  invisible(beta)
}

ocCurves.xbar <- function(object, n, c = seq(0, 5, length=101), 
                          nsigmas = object$nsigmas, 
                          lty = rep(1,length(n)),
                          lwd = rep(2,length(n)),
                          col = blues.colors(length(n)),
                          identify=FALSE, 
                          restore.par=TRUE)
{
# Draw the operating-characteristic curves for the xbar-chart with nsigmas
# limits. The values on the vertical axis give the probability of not detecting
# a shift of c*sigma in the mean on the first sample following the shift.

  if (!(object$type == "xbar"))
     stop("not a 'qcc' object of type \"xbar\".")

  size <- unique(object$sizes)
  if (length(size) > 1)
     stop("Operating characteristic curves available only for equal sample sizes!")
  if (missing(n))
     n <- sort(unique(c(size, c(1,5,10,15,20))))
  if (is.null(nsigmas))
     nsigmas <- qnorm(1 - (1 - object$confidence.level) / 2)

  beta <- matrix(as.double(NA), length(c), length(n))
  for(i in 1:length(n))
     beta[,i] <- pnorm(nsigmas-c*sqrt(n[i])) - pnorm(-nsigmas-c*sqrt(n[i]))
  colnames(beta) <- paste("n=",n,sep="")
  rownames(beta) <- c
  names(dimnames(beta)) <- c("shift (std.dev)", "sample size")

  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  cex.labels <- par("cex")*qcc.options("cex")
  cex.stats <- par("cex")*qcc.options("cex.stats")
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = c(4.1,4.1,1.1,2.1),
      oma = c(0, 0, 1.5*cex.labels, 0))

  plot(c, beta[,1], type="n",
       ylim = c(0,1), xlim = c(0,max(c)),
       xlab = "Process shift (std.dev)",
       ylab = "Prob. type II error ")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  box()
  mtext(paste("OC curves for", object$type, "chart"), 
        side = 3, outer = TRUE, line = 0, adj = 0, at = par("plt")[1],
        font = par("font.main"), 
        cex  = par("cex")*qcc.options("cex"), 
        col  = par("col.main"))
  matlines(c, beta, lty = lty, lwd = lwd, col = col)
  
  if(identify)
  { 
    cs <- rep(c,length(n))
    betas <- as.vector(beta)
    labels <- paste("c=", formatC(cs, 2, flag="-"), 
                    ": beta=", formatC(betas, 4, flag="-"), 
                    ", ARL=", formatC(1/(1-betas), 2, flag="-"), sep="")
    i <- identify(cs, betas, labels, pos=4, offset=0.2,
                  cex = cex.stats)
    apply(as.matrix(labels[i$ind]), 1, cat, "\n")
  } else
  { 
    legend("topright", inset = 0.02, 
           legend = paste("n =", n), cex = cex.stats,
           bg = qcc.options("bg.figure"),
           lty = lty, lwd = lwd, col = col)
  }
  
  invisible(beta)
}

ocCurves.p <- function(object, nsigmas = object$nsigmas, 
                        lty = 1, lwd = 2, col = blues.colors(1),
                        identify = FALSE, restore.par = TRUE)
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

  if (object$type=="p") 
     { UCL <- min(floor(size*limits[,2]), size)
       LCL <- max(floor(size*limits[,1]), 0) }
  else
     { UCL <- min(floor(limits[,2]), size)
       LCL <- max(floor(limits[,1]), 0) }
  beta <- pbinom(UCL, size, p) - pbinom(LCL-1, size, p)
  names(beta) <- p

  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  cex.labels <- par("cex")*qcc.options("cex")
  cex.stats <- par("cex")*qcc.options("cex.stats")
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = c(4.1,4.1,1.1,2.1),
      oma = c(0, 0, 1.5*cex.labels, 0))
  
  plot(p, beta, type = "n", 
       ylim = c(0,1), xlim = c(0,1),
       xlab = expression(p), 
       ylab = "Prob. type II error ")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  box()
  mtext(paste("OC curves for", object$type, "chart"), 
        side = 3, outer = TRUE, line = 0, adj = 0, at = par("plt")[1],
        font = par("font.main"), 
        cex  = par("cex")*qcc.options("cex"), 
        col  = par("col.main"))
  lines(p, beta, lty = lty, lwd = lwd, col = col)
  lines(rep(p[which.max(beta)], 2), c(0, max(beta)), lty = 2)

  warning("Some computed values for the type II error have been rounded due to the discreteness of the binomial distribution. Thus, some ARL values might be meaningless.")

  if (identify)
     { labels <- paste("p=", formatC(p, 2, flag="-"), 
                       ": beta=", formatC(beta, 4, flag="-"), 
                       ", ARL=", formatC(1/(1-beta), 2, flag="-"), sep="")
       i <- identify(p, beta, labels, pos=4, offset=0.2,
                     cex = cex.stats)
       apply(as.matrix(labels[i$ind]), 1, cat, "\n")
  }
  
  invisible(beta)  
}

ocCurves.c <- function(object, nsigmas = object$nsigmas, 
                        lty = 1, lwd = 2, col = blues.colors(1),
                        identify = FALSE, restore.par = TRUE)
{
  type <- object$type
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
  if (object$type=="c") 
     { max.lambda <- ceiling(CL+10*std.dev)
       UCL <- floor(limits[1,2])
       LCL <- floor(limits[1,1])
     }
  else
     { max.lambda <- ceiling(CL*size+10*std.dev)[1]
       UCL <- floor(size*limits[1,2])
       LCL <- floor(size*limits[1,1])
     }
  lambda <- seq(0, max.lambda)
  beta <- ppois(UCL, lambda) - ppois(LCL-1, lambda)
  names(beta) <- lambda

  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  cex.labels <- par("cex")*qcc.options("cex")
  cex.stats <- par("cex")*qcc.options("cex.stats")
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = c(4.1,4.1,1.1,2.1),
      oma = c(0, 0, 1.5*cex.labels, 0))

  plot(lambda, beta, type = "n", 
       ylim = c(0,1), xlim = range(lambda),
       xlab = "Mean", 
       ylab = "Prob. type II error ")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  box()
  mtext(paste("OC curves for", object$type, "chart"), 
        side = 3, outer = TRUE, line = 0, adj = 0, at = par("plt")[1],
        font = par("font.main"), 
        cex  = qcc.options("cex"), 
        col  = par("col.main"))
  lines(lambda, beta, lty = lty, lwd = lwd, col = col)
  lines(rep(lambda[which.max(beta)], 2), c(0, max(beta)), lty = 2)

  warning("Some computed values for the type II error have been rounded due to the discreteness of the Poisson distribution. Thus, some ARL values might be meaningless.")

  if (identify)
     { labels <- paste("lambda=", formatC(lambda, 0, flag="-"), 
                       ": beta=", formatC(beta, 4, flag="-"), 
                       ", ARL=", formatC(1/(1-beta), 2, flag="-"), sep="")
       i <- identify(lambda, beta, labels, pos=4, offset=0.2,
                     cex = cex.stats)
       apply(as.matrix(labels[i$ind]), 1, cat, "\n")
  }
  
  invisible(beta)
}

ocCurves.R <- function(object, n, c = seq(1, 6, length=101), 
                        nsigmas = object$nsigmas, 
                        lty = rep(1,length(n)),
                        lwd = rep(2,length(n)),
                        col = blues.colors(length(n)),
                        identify=FALSE, 
                        restore.par=TRUE)
{
# Draw the operating-characteristic curves for the R-chart with nsigmas
# limits. The values on the vertical axis give the probability of not detecting
# a change from sigma to c*sigma on the first sample following the change.

  if (!(object$type=="R"))
     stop("not a `qcc' object of type \"R\".")

  size <- unique(object$sizes)
  if (length(size) > 1)
     stop("Operating characteristic curves available only for equal sample sizes!")
  if (missing(n))
     n <- sort(unique(c(size, c(2,5,10,15,20))))
  if (is.null(nsigmas))
    { tail.prob <- (1 - object$confidence.level) / 2
      beta.fun1 <- function(c, n, p)
      {
        lcl <- qtukey(p, n, Inf)
        ucl <- qtukey(p, n, Inf, lower.tail = FALSE)
        ptukey(ucl / c, n, Inf) - ptukey(lcl / c, n, Inf)
      }
      beta <- outer(c, n, beta.fun1, tail.prob)
    }
  else
    { exp.R.unscaled <- qcc.options("exp.R.unscaled")
      se.R.unscaled <- qcc.options("se.R.unscaled")
      Rtab <- min(length(exp.R.unscaled), length(se.R.unscaled))
      if (any(n > Rtab))
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
      beta <- outer(c, n, beta.fun2, nsigmas)
    }
  colnames(beta) <- paste("n=",n,sep="")
  rownames(beta) <- c
  names(dimnames(beta)) <- c("scale multiplier", "sample size")

  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  cex.labels <- par("cex")*qcc.options("cex")
  cex.stats <- par("cex")*qcc.options("cex.stats")
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = c(4.1,4.1,1.1,2.1),
      oma = c(0, 0, 1.5*cex.labels, 0))

  plot(c, beta[,1], type="n",
       ylim = c(0,1), xlim = c(1,max(c)),
       xlab = "Process scale multiplier",
       ylab = "Prob. type II error ")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
       col = qcc.options("bg.figure"))
  box()
  mtext(paste("OC curves for", object$type, "chart"), 
        side = 3, outer = TRUE, line = 0, adj = 0, at = par("plt")[1],
        font = par("font.main"), 
        cex  = qcc.options("cex"), 
        col  = par("col.main"))
  matlines(c, beta, lty = lty, lwd = lwd, col = col)

  if(identify)
  { 
    cs <- rep(c,length(n))
    betas <- as.vector(beta)
    labels <- paste("c=", formatC(cs, 2, flag="-"),
                    ": beta=", formatC(betas, 4, flag="-"),
                    ", ARL=", formatC(1/(1-betas), 2, flag="-"), sep="")
    i <- identify(cs, betas, labels, pos=4, offset=0.2,
                  cex = cex.stats)
    apply(as.matrix(labels[i$ind]), 1, cat, "\n")
  } else
  { 
    legend("topright", inset = 0.02, 
           legend = paste("n =", n), cex = cex.stats,
           bg = qcc.options("bg.figure"),
           lty = lty, lwd = lwd, col = col)
  }
  
  invisible(beta)
}

ocCurves.S <- function(object, n, c = seq(1, 6, length=101), 
                        nsigmas = object$nsigmas, 
                        lty = rep(1,length(n)),
                        lwd = rep(2,length(n)),
                        col = blues.colors(length(n)),
                        identify = FALSE, 
                        restore.par = TRUE)
{
# Draw the operating-characteristic curves for the S-chart with nsigmas
# limits. The values on the vertical axis give the probability of not detecting
# a change from sigma to c*sigma on the first sample following the change.

  if (!(object$type=="S"))
     stop("not a `qcc' object of type \"S\".")

  size <- unique(object$sizes)
  if (length(size) > 1)
     stop("Operating characteristic curves available only for equal sample sizes!")
  if (missing(n))
     n <- sort(unique(c(size, c(2,5,10,15,20))))
  if (is.null(nsigmas))
    { tail.prob <- (1 - object$confidence.level) / 2
      beta.fun1 <- function(c, n, p)
      {
        ucl <- sqrt(qchisq(1 - p, n - 1) / (n - 1))
        lcl <- sqrt(qchisq(p, n - 1) / (n - 1))
        pchisq((n - 1) * (ucl / c)^2, n - 1) - pchisq((n - 1)* (lcl / c)^2, n - 1)
      }
      beta <- outer(c, n, beta.fun1, tail.prob)
    }
  else
    { c4 <- qcc.c4
      beta.fun2 <- function(c, n)
      {
        center <- c4(n)
        tol <- sqrt(1 - c4(n)^2)
        lcl <- pmax(0, center - nsigmas * tol)
        ucl <- center + nsigmas * tol
        pchisq((n - 1) * (ucl / c)^2, n - 1) - pchisq((n - 1) * (lcl / c)^2, n - 1)
      }
      beta <- outer(c, n, beta.fun2)
    }
  colnames(beta) <- paste("n=",n,sep="")
  rownames(beta) <- c
  names(dimnames(beta)) <- c("scale multiplier", "sample size")

  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  cex.labels <- par("cex")*qcc.options("cex")
  cex.stats <- par("cex")*qcc.options("cex.stats")
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = c(4.1,4.1,1.1,2.1),
      oma = c(0, 0, 1.5*cex.labels, 0))
  
  plot(c, beta[,1], type="n",
       ylim = c(0,1), xlim = c(1,max(c)),
       xlab = "Process scale multiplier",
       ylab = "Prob. type II error ")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
       col = qcc.options("bg.figure"))
  box()
  mtext(paste("OC curves for", object$type, "Chart"), 
        side = 3, outer = TRUE, line = 0, adj = 0, at = par("plt")[1],
        font = par("font.main"), 
        cex  = qcc.options("cex"), 
        col  = par("col.main"))
  matlines(c, beta, lty = lty, lwd = lwd, col = col)

  if (identify)
     { cs <- rep(c,length(n))
       betas <- as.vector(beta)
       labels <- paste("c=", formatC(cs, 2, flag="-"),
                       ": beta=", formatC(betas, 4, flag="-"),
                       ", ARL=", formatC(1/(1-betas), 2, flag="-"), sep="")
       i <- identify(cs, betas, labels, pos=4, offset=0.2,
                     cex = cex.stats)
       apply(as.matrix(labels[i$ind]), 1, cat, "\n")
     }
  else
     { legend("topright", inset = 0.02, 
              legend = paste("n =", n), cex = cex.stats,
              bg = qcc.options("bg.figure"),
              lty = lty, lwd = lwd, col = col)
     }
  
  invisible(beta)
}
