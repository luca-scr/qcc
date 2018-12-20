#-------------------------------------------------------------------#
#                                                                   #
#                      CUSUM CHART                                  #
#                                                                   #
#-------------------------------------------------------------------#

cusum <- function(data, sizes, center, std.dev, head.start = 0, 
                  decision.interval = 5, se.shift = 1,
                  data.name, labels, newdata, newsizes, newlabels, 
                  plot = TRUE, ...)
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

  if (decision.interval <= 0)
      stop("decision.interval must be positive")

  if (head.start < 0 || head.start >= decision.interval)
      stop("head.start must be non-negative and less than decision.interval")
  
  # used for computing statistics and std.dev
  type <- if(any(sizes==1)) "xbar.one" else "xbar"
  if(ncol(data) == 1 & any(sizes > 1) & missing(std.dev))
     stop("sizes larger than 1 but data appears to be single samples. In this case you must provide also the std.dev")
  
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

  stopifnot(length(labels) == length(statistics))
  names(statistics) <-  rownames(data) <-  labels
  names(dimnames(data)) <- list("Group", "Samples")

  object <- list(call = call, type = "cusum", 
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
      stopifnot(length(newlabels) == length(newstats))
      names(newstats) <- newlabels
      object$newstats <- newstats
      object$newdata  <- newdata
      object$newsizes <- newsizes
      object$newdata.name <- newdata.name
      statistics <- c(statistics, newstats)
      sizes <- c(sizes, newsizes)
    }

  n <- length(statistics)
  z <- (statistics - center)/(std.dev/sqrt(sizes))
  ldb <-  - decision.interval
  udb <- decision.interval
  #
  z.f <- z - se.shift/2
  cusum.pos <- rep(NA, n)
  cusum.pos[1] <- max(0, head.start + z.f[1])
  for (i in 2:n)
      cusum.pos[i] <- max(0, cusum.pos[i-1]+z.f[i])
  #
  z.f <- z + se.shift/2
  cusum.neg <- rep(NA, n)
  cusum.neg[1] <- max(0, head.start - z.f[1])
  for (i in 2:n)
      cusum.neg[i] <- max(0, cusum.neg[i-1]-z.f[i])
  cusum.neg <- -cusum.neg
  violations <- list(lower = which(cusum.neg < ldb),
                     upper = which(cusum.pos > udb))
  
  object$type <- "cusum"
  object$pos <- cusum.pos 
  object$neg <- cusum.neg
  object$head.start <- head.start
  object$decision.interval <- decision.interval
  object$se.shift <- se.shift
  object$violations <- violations

  class(object) <- "cusum.qcc"
  if(plot) plot(object, ...) 
  return(object)
}

print.cusum.qcc <- function(x, digits =  getOption("digits"), ...)
{
  object <- x   # Argh.  Really want to use 'object' anyway
  # cat("\nCall:\n",deparse(object$call),"\n\n",sep="")
  cat(cli::rule(left = crayon::bold("Cusum Chart"), 
                width = min(getOption("width"),50)), "\n\n")
  
  data.name <- object$data.name
  # type <- object$type
  statistics <- object$statistics
  # cat(paste(type, "chart for", data.name, "\n"))
  # cat("\nSummary of group statistics:\n")
  # print(summary(statistics), ...)
  
  cat("Data (phase I)             =", data.name, "\n")
  cat("Number of groups           =", length(statistics), "\n")
  
  sizes <- object$sizes
  if(length(unique(sizes))==1)
     sizes <- sizes[1]
  if(length(sizes) == 1)
  {
    cat("Group sample size          =", signif(sizes), "\n")
  } else 
  {
    cat("Group sample sizes         =")
    tab <- table(sizes)
    print(matrix(c(as.numeric(names(tab)), tab), 
                 ncol = length(tab), byrow = TRUE, 
                 dimnames = list(c("  sizes", "  counts"),
                                 character(length(tab)))),
          digits = digits, ...)
  }

  center <- object$center
  cat("Center of group statistics =", signif(center, digits = digits), "\n")

  sd <- object$std.dev
  cat("Standard deviation         =", signif(sd, digits = digits), "\n")

  newdata.name <- object$newdata.name
  newstats <- object$newstats
  if(!is.null(newstats)) 
  { 
    # cat(paste("\nSummary of group statistics in ", 
    #           newdata.name, ":\n", sep = ""))
    # print(summary(newstats), digits = digits, ...)
    cat("\nNew data (phase II)        =", newdata.name, "\n")
    cat("Number of groups           =", length(newstats), "\n")
    newsizes <- object$newsizes
    if (length(unique(newsizes)) == 1)
      newsizes <- newsizes[1]
    if(length(newsizes) == 1)
    {
      cat("Group sample size          =", signif(newsizes), "\n")
    } else 
    { 
      cat("Group sample sizes         =")
      new.tab <- table(newsizes)
      print(matrix(c(as.numeric(names(new.tab)), new.tab),
                   ncol = length(new.tab), byrow = TRUE, 
                   dimnames = list(c("  sizes", "  counts"),
                                   character(length(new.tab)))),
            digits = digits, ...)
    }
  }
  
  cat("\n")
  if(object$head.start > 0)
  {
    cat("Head start (std.err.)      =",
        signif(object$head.start, digits = digits), "\n")
  }
  cat("Decision interval (std.err.) =", 
      signif(object$decision.interval, digits = digits), "\n")
  cat("Shift detection   (std.err.) =", 
      signif(object$se.shift, digits = digits), "\n")

  invisible()
}

summary.cusum.qcc <- function(object, ...) print.cusum.qcc(object, ...)

plot.cusum.qcc <- function(x, 
                           add.stats = qcc.options("add.stats"), 
                           chart.all = qcc.options("chart.all"), 
                           fill = qcc.options("fill"),
                           label.bounds = c("LDB", "UDB"), 
                           title, xlab, ylab, ylim, axes.las = 0,
                           digits = getOption("digits"),
                           restore.par = TRUE, ...) 
{
  object <- x  # Argh.  Really want to use 'object' anyway
  if ((missing(object)) | (!inherits(object, "cusum.qcc")))
     stop("an object of class `cusum.qcc' is required")

  # collect info from object
  type <- object$type
  data.name <- object$data.name
  center <- object$center
  std.dev <- object$std.dev
  stats <- object$statistics
  ldb <- -object$decision.interval
  udb <- object$decision.interval
  newstats <- object$newstats
  newdata.name <- object$newdata.name
  violations <- object$violations
  cusum.pos <- object$pos
  cusum.neg <- object$neg

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
         main.title <- paste(type, "chart for", data.name)
      else if(chart.all)
                main.title <- paste(type, "chart for", data.name, 
                                  "and", newdata.name)
           else main.title <- paste(type, "Chart for", newdata.name) }
  else main.title <- paste(title)
  cex.labels <- par("cex")*qcc.options("cex")
  cex.stats <- par("cex")*qcc.options("cex.stats")

  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      # mgp = c(2.1, 0.8, 0),
      mar = c(4.1,4.1,1.1,2.1),
      oma = if(add.stats) c(3.5*cex.stats, 0, 1.5*cex.labels, 0) 
            else          c(0, 0, 1.5*cex.labels, 0))

  plot(indices, statistics, type="n",
       ylim = if(!missing(ylim)) ylim 
              else range(cusum.pos, cusum.neg, ldb, udb),
       ylab = "", xlab = if(missing(xlab)) "Group" else xlab, 
       axes = FALSE)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  axis(1, at = indices, las = axes.las,
       labels = if(is.null(names(statistics))) 
                   as.character(indices) else names(statistics),
       cex.axis = par("cex.axis")*0.9)
  axis(2, las = axes.las, cex.axis = par("cex.axis")*0.9)
  box()
  mtext(main.title, side = 3, outer = TRUE,
        line = 0, adj = 0, at = par("plt")[1],
        font = par("font.main"), 
        cex  = par("cex")*qcc.options("cex"), 
        col  = par("col.main"))

  mtext(if(missing(ylab)) "Cumulative Sum" else ylab, line=3, side=2)
  lab <- "Above target"
  if (add.stats && object$head.start > 0)
      lab <- paste(lab, " (start = ", object$head.start, ")", sep = "")
  mtext(lab, srt=90, line=2, side=2, at=0+par("usr")[4]/2,
        cex = par("cex")*0.8)
  lab <- "Below target"
  if (add.stats && object$head.start > 0)
      lab <- paste(lab, " (start = ", - object$head.start, ")", sep = "")
  mtext(lab, srt=90, line=2, side=2, at=0+par("usr")[3]/2,
        cex = par("cex")*0.8)
  mtext(label.bounds, side = 4, at = c(ldb, udb), las = 1, line = 0.1, 
        col = gray(0.3), cex = par("cex")*qcc.options("cex.stats"))

  # draw decision area/lines
  if(fill)
  { 
    polygon(par("usr")[c(1,2,2,1)], c(ldb, ldb, udb, udb), 
            border = FALSE, 
            col = adjustcolor(qcc.options("zones")$fill, alpha.f = 0.1)) 
  } else
  {
    abline(h = c(ldb, udb), 
           lty = qcc.options("zones")$lty[1], 
           col = qcc.options("zones")$col[1])
  }
  
  # draw lines & points
  abline(h = 0, col = qcc.options("zones")$col[1])
  lines(indices, cusum.pos[indices], type = "b", pch=NA)
  lines(indices, cusum.neg[indices], type = "b", pch=NA)
  #
  col <- rep(palette()[1], length(cusum.pos))
  pch <- rep(20, length(cusum.pos))
  if(!is.null(violations$upper))
  { 
    col[violations$upper] <- qcc.options("rules")$col[1] 
    pch[violations$upper] <- qcc.options("rules")$pch[1]  
  }
  points(indices, cusum.pos[indices], col = col[indices], pch = pch[indices])
  #
  col <- rep(palette()[1], length(cusum.neg))
  pch <- rep(20, length(cusum.neg))
  if(!is.null(violations$lower))
  { 
    col[violations$lower] <- qcc.options("rules")$col[1] 
    pch[violations$lower] <- qcc.options("rules")$pch[1]  
  }
  points(indices, cusum.neg[indices], col = col[indices], pch = pch[indices])

  if(chart.all & !is.null(newstats))
  { 
    len.obj.stats <- length(stats)
    len.new.stats <- length(newstats)
    abline(v = len.obj.stats + 0.5, lty = 3)
    mtext("Calibration data", cex = par("cex")*0.8, 
          at = len.obj.stats/2, line = 0, adj = 0.5)
    mtext("New data", cex = par("cex")*0.8, 
          at = len.obj.stats + len.new.stats/2, line = 0, adj = 0.5)
  }

  if(add.stats) 
  { 
    at <- c(0.15,0.55) 
    # write info at bottom
    mtext(paste("Number of groups = ", length(statistics), sep = ""),
          side = 1, outer = TRUE, line = 0, adj = 0, at = at[1],
          font = qcc.options("font.stats"),
          cex = par("cex")*qcc.options("cex.stats"))
    if(length(center) == 1)
    { mtext(paste("Center = ", signif(center[1], digits), sep = ""),
            side = 1, outer= TRUE, line = 1*cex.stats, adj = 0, at = at[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    } else 
    { mtext("Center is variable", 
            side = 1, outer = TRUE, line = 1*cex.stats, adj = 0, at = at[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    }
    mtext(paste("StdDev = ", signif(std.dev, digits), sep = ""),
          side = 1, outer = TRUE, line = 2*cex.stats, adj = 0, at = at[1],
          font = qcc.options("font.stats"),
          cex = par("cex")*qcc.options("cex.stats"))
    mtext(paste("Decision interval (std. err.) =", 
                signif(object$decision.interval, digits = digits)),
          side = 1, outer = TRUE, line = 0, adj = 0, at = at[2],
          font = qcc.options("font.stats"),
          cex = par("cex")*qcc.options("cex.stats"))
    mtext(paste("Shift detection (std. err.) =", 
                signif(object$se.shift, digits = digits)), 
          side = 1, outer = TRUE, line = 1*cex.stats, adj = 0, at = at[2], 
          font = qcc.options("font.stats"),
          cex = par("cex")*qcc.options("cex.stats"))
    mtext(paste("No. of points beyond boundaries =", 
                sum(sapply(violations, length))),
          side = 1, outer = TRUE, line = 2*cex.stats, adj = 0, at = at[2],
          font = qcc.options("font.stats"),
          cex = par("cex")*qcc.options("cex.stats"))
  }

  invisible()
}
