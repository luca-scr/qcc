#-------------------------------------------------------------------#
#                                                                   #
#                       EWMA CHART                                  #
#                                                                   #
#-------------------------------------------------------------------#

ewmaSmooth <- function(x, y, lambda = 0.20, start, ...)
{
#
# Exponential-Weighted Moving Average 
# 
# Return smooth values based on 
# 
# z_t = lambda*y_t + (1-lambda)*z_t-1      
# 
# where 0<= lambda <=1 is the parameter which controls the weights applied 
# to the data, and start is the starting value.
# Returns a list with elements:
# x = ordered x-values
# y = smoothed fitted values of y
# 
  if (length(y)!=length(x))
     stop("x and y must have the same length!")
  if (abs(lambda)>1)
     stop("lambda parameter must be between 0 and 1")
  ord <- order(x) 
  x <- x[ord]
  y <- y[ord]
  n <- length(y)
  if (missing(start)) start <- y[1]
  z <- c(start, y)
  for (i in 2:(n + 1))
    z[i] <- lambda * z[i] + (1 - lambda) * z[i - 1]
  list(x=x, y=z[-1], lambda=lambda, start=start)
}


ewma <- function(data, sizes, center, std.dev, lambda = 0.2, nsigmas = 3, data.name, labels, newdata, newsizes, newlabels, plot = TRUE, ...)
{

  call <- match.call()
  if (missing(data))
     stop("'data' argument is not specified")

  if(missing(data.name)) 
     data.name <- deparse(substitute(data))
  data <- data.matrix(data)

  if(missing(sizes)) 
    { sizes <- apply(data, 1, function(x) sum(!is.na(x)))  }
  else
    { if(length(sizes)==1)
         sizes <- rep(sizes, nrow(data))
      else if(length(sizes) != nrow(data))
              stop("sizes length doesn't match with data") }
  # used for computing statistics and std.dev
  type <- if(any(sizes==1)) "xbar.one" else "xbar"

  if(missing(labels))
    { if(is.null(rownames(data))) labels <- 1:nrow(data)
      else                        labels <- rownames(data) }

  stats <- paste("stats.", type, sep = "")
  if(!exists(stats, mode="function"))
     stop(paste("function", stats, "is not defined"))
  stats <- do.call(stats, list(data, sizes))
  statistics <- stats$statistics
  if(missing(center)) center <- stats$center

  sd <- paste("sd.", type, sep = "")
  if(!exists(sd, mode="function"))
     stop(paste("function", sd, "is not defined!"))
  if(missing(std.dev)) 
    { std.dev <- switch(type, 
                        "xbar" = { if(any(sizes > 25)) "RMSDF"
                                   else                "UWAVE-R" },
                         NULL)
      std.dev <- do.call(sd, list(data, sizes, std.dev)) }
  else 
     { if (is.character(std.dev))
          { std.dev <- do.call(sd, list(data, sizes, std.dev)) }
       else
          { if (!is.numeric(std.dev))
               stop("if provided the argument 'std.dev' must be a method available or a numerical value. See help(qcc).")  }
     }

  names(statistics) <-  rownames(data) <-  labels
  names(dimnames(data)) <- list("Group", "Samples")
  object <- list(call = call, type = "ewma", 
                 data.name = data.name, data = data, 
                 statistics = statistics, sizes = sizes, 
                 center = center, std.dev = std.dev)

  # check for new data provided and update object
  if(!missing(newdata))
    { newdata.name <- deparse(substitute(newdata))
      newdata <- data.matrix(newdata)
      if(missing(newsizes))
        { newsizes <- apply(newdata, 1, function(x) sum(!is.na(x))) }
      else
        { if(length(newsizes)==1)
            newsizes <- rep(newsizes, nrow(newdata))
          else if(length(newsizes) != nrow(newdata))
                  stop("newsizes length doesn't match with newdata") }
      stats <- paste("stats.", type, sep = "")
      if(!exists(stats, mode="function"))
         stop(paste("function", stats, "is not defined"))
      newstats <- do.call(stats, list(newdata, newsizes))$statistics
      if(missing(newlabels))
        { if(is.null(rownames(newdata)))
            { start <- length(statistics)
              newlabels <- seq(start+1, start+length(newstats)) }
          else
            { newlabels <- rownames(newdata) }
        }
      names(newstats) <- newlabels
      object$newstats <- newstats
      object$newdata  <- newdata
      object$newsizes <- newsizes
      object$newdata.name <- newdata.name
      statistics <- c(statistics, newstats)
      sizes <- c(sizes, newsizes)
    }

  n <- length(statistics)
  indices <- 1:length(statistics)
  ewma <- ewmaSmooth(indices, statistics, lambda=lambda, start=center)
  sigma2 <- std.dev^2/sizes * 
            ((lambda/(2-lambda))*(1-(1-lambda)^(2*(1:n))))
  ucl <- center + nsigmas*sqrt(sigma2)
  lcl <- center - nsigmas*sqrt(sigma2)

  object$x <- ewma$x
  y <- as.vector(ewma$y)
  names(y) <- c(names(object$statistics), names(object$newstats))
  object$y <- y
  object$sigma <- sqrt(sigma2)
  object$lambda <- lambda
  object$nsigmas <- nsigmas
  limits <- cbind(lcl,ucl)
  colnames(limits) <- c("LCL", "UCL")
  object$limits <- limits
  object$violations <- which(y < lcl | y > ucl)

  class(object) <- "ewma.qcc"
  if(plot) plot(object, ...) 
  return(object)
}


print.ewma.qcc <- function(x, ...) str(x, 1)

summary.ewma.qcc <- function(object, digits =  getOption("digits"), ...)
{
  # object <- x   # Argh.  Really want to use 'object' anyway
  cat("\nCall:\n",deparse(object$call),"\n\n",sep="")
  data.name <- object$data.name
  type <- object$type
  cat(paste(type, "chart for", data.name, "\n"))
  statistics <- object$statistics
  cat("\nSummary of group statistics:\n")
  print(summary(statistics), digits = digits, ...)
  sizes <- object$sizes
  if(length(unique(sizes))==1)
     sizes <- sizes[1]
  if(length(sizes) == 1)
     cat("\nGroup sample size: ", format(sizes))
  else {
         cat("\nSummary of group sample sizes: ")
         tab <- table(sizes)
         print(matrix(c(as.numeric(names(tab)), tab), 
                      ncol = length(tab), byrow = TRUE, 
                      dimnames = list(c("  sizes", "  counts"),
                      character(length(tab)))), 
               digits = digits, ...)
        }
  cat("\nNumber of groups: ", length(statistics))
  center <- object$center
  cat("\nCenter of group statistics: ", format(center, digits = digits))
  sd <- object$std.dev
  cat("\nStandard deviation: ", format(sd, digits = digits), "\n")

  newdata.name <- object$newdata.name
  newstats <- object$newstats
  if(!is.null(newstats)) 
    { cat(paste("\nSummary of group statistics in ", 
                newdata.name, ":\n", sep = ""))
      print(summary(newstats), digits = digits, ...)
      newsizes <- object$newsizes
      if (length(unique(newsizes)) == 1)
         newsizes <- newsizes[1]
      if(length(newsizes) == 1)
         cat("\nGroup sample size: ", format(newsizes))
      else 
         { cat("\nSummary of group sample sizes:\n")
           new.tab <- table(newsizes)
           print(matrix(c(as.numeric(names(new.tab)), new.tab),
                        ncol = length(new.tab), byrow = TRUE, 
                        dimnames = list(c("  sizes", "  counts"),
                                        character(length(new.tab)))), 
                 digits = digits, ...)
         }
      cat("\nNumber of groups: ", length(newstats), "\n")
     }

  cat("\nSmoothing parameter:", signif(object$lambda, digits = digits), "\n")

  limits <- object$limits
  if (!is.null(limits)) 
     { cat("Control limits:\n")
       .printShortMatrix(limits, digits = digits, ...) }
       
  invisible()
}

plot.ewma.qcc <- function(x, add.stats = TRUE, chart.all = TRUE, 
                          label.limits = c("LCL", "UCL"),
                          title, xlab, ylab, ylim, axes.las = 0,
                          digits =  getOption("digits"),
                          restore.par = TRUE, ...) 
{
  object <- x  # Argh.  Really want to use 'object' anyway
  if((missing(object)) | (!inherits(object, "ewma.qcc")))
     stop("an object of class `ewma.qcc' is required")

  # collect info from object
  type <- object$type
  data.name <- object$data.name
  center <- object$center
  std.dev <- object$std.dev
  stats <- object$statistics
  limits <- object$limits
  newstats <- object$newstats
  newdata.name <- object$newdata.name
  violations <- object$violations
  nviolations <- length(violations)

  if(chart.all) 
    { statistics <- c(stats, newstats) 
      indices <- 1:length(statistics) }
  else
    { if(is.null(newstats))
        { statistics <- stats 
          indices <- 1:length(statistics) }
      else
         { statistics <- newstats
           indices <- seq(length(stats)+1, length(stats)+length(newstats)) }
    } 
     
  if(missing(title))
    { if(is.null(newstats))
         main.title <- paste("EWMA Chart\nfor", data.name)
      else if(chart.all)
                main.title <- paste("EWMA Chart\nfor", data.name, 
                                  "and", newdata.name)
           else main.title <- paste("EWMA Chart\nfor", newdata.name) }
  else main.title <- paste(title)

  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  mar <- pmax(oldpar$mar, c(4.1,4.1,3.1,2.1))
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex*qcc.options("cex"),
      mar = if(add.stats) pmax(mar, c(7.6,0,0,0)) else mar)

  plot(indices, statistics, type="n",
       ylim = if(!missing(ylim)) ylim 
              else range(statistics, limits),
       ylab = ifelse(missing(ylab), "Group Summary Statistics", ylab),
       xlab = ifelse(missing(xlab), "Group", xlab),
       axes = FALSE)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  axis(1, at = indices, las = axes.las,
       labels = if(is.null(names(statistics))) 
                   as.character(indices) else names(statistics))
  axis(2, las = axes.las)
  box()
  top.line <- par("mar")[3] - length(capture.output(cat(main.title)))
  top.line <- top.line - if(chart.all & (!is.null(newstats))) 0.1 else 0.5
  mtext(main.title, side = 3, line = top.line,
        font = par("font.main"), 
        cex  = qcc.options("cex"), 
        col  = par("col.main"))

  abline(h = center, lty=1)
  lines(indices, limits[indices,1], lty=2)
  lines(indices, limits[indices,2], lty=2)

  points(indices, statistics, pch = 3, cex = 0.8)
  lines(indices, object$y[indices], type = "o", pch=20)
  
  mtext(label.limits, side = 4, at = limits[nrow(limits),], 
        las = 1, line = 0.1, col = gray(0.3), cex = par("cex"))

  if(nviolations > 0)
    { if(is.null(qcc.options("beyond.limits")))
         stop(".qcc.options$beyond.limits undefined. See help(qcc.options).")
      points(violations, object$y[violations], 
             col = qcc.options("beyond.limits")$col, 
             pch = qcc.options("beyond.limits")$pch)
    }
              
  if(chart.all & !is.null(newstats))
    { len.obj.stats <- length(stats)
      len.new.stats <- length(newstats)
      abline(v = len.obj.stats + 0.5, lty = 3)
      mtext(# paste("Calibration Data in", data.name), 
            "Calibration data", cex = par("cex")*0.8,
            at = len.obj.stats/2, line = 0, adj = 0.5)
      mtext(# paste("New Data in", object$newdata.name), 
            "New data", cex = par("cex")*0.8, 
            at = len.obj.stats + len.new.stats/2, line = 0, adj = 0.5)
     }

  if(add.stats) 
    { # computes the x margins of the figure region
      plt <- par()$plt; usr <- par()$usr
      px <- diff(usr[1:2])/diff(plt[1:2])
      xfig <- c(usr[1]-px*plt[1], usr[2]+px*(1-plt[2]))
      at.col <- xfig[1] + diff(xfig[1:2])*c(0.15, 0.6)
      top.line <- 4.5
      # write info at bottom
      mtext(paste("Number of groups = ", length(statistics), sep = ""),
           side = 1, line = top.line, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      if(length(center) == 1)
         { mtext(paste("Center = ", signif(center[1], digits = digits), sep = ""),
                 side = 1, line = top.line+1, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
         }
       else 
         { mtext("Center is variable", 
                 side = 1, line = top.line+1, adj = 0, at = at.col[1],
                 font = qcc.options("font.stats"),
                 cex = par("cex")*qcc.options("cex.stats"))
         }
      mtext(paste("StdDev = ", signif(std.dev, digits = digits), sep = ""),
            side = 1, line = top.line+2, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Smoothing parameter = ", signif(object$lambda, digits = digits)),
            side = 1, line = top.line, adj = 0, at = at.col[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Control limits at ", object$nsigmas, "*sigma", sep=""),  
            side = 1, line = top.line+1, adj = 0, at = at.col[2], 
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("No. of points beyond limits =", nviolations),
            side = 1, line = top.line+2, adj = 0, at = at.col[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      }

  invisible()
}
