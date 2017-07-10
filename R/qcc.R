#-----------------------------------------------------------------------------#
#                                                                             #
#                     QUALITY CONTROL CHARTS IN R                             #
#                                                                             #
#  An R package for statistical in-line quality control.                      #
#                                                                             #
#  Written by: Luca Scrucca                                                   #
#              Department of Statistics                                       #
#              University of Perugia, ITALY                                   #
#              luca@stat.unipg.it                                             #
#                                                                             #
#-----------------------------------------------------------------------------#

#
#  Main function to create a 'qcc' object
#

qcc <- function(data, type = c("xbar", "R", "S", "xbar.one", "p", "np", "c", "u", "g"), sizes, center, std.dev, limits, data.name, labels, newdata, newsizes, newdata.name, newlabels, nsigmas = 3, confidence.level, rules = shewhart.rules, plot = TRUE, ...)
{
  call <- match.call()
  
  if (missing(data))
     stop("'data' argument is not specified")
  
  if(identical(type, eval(formals(qcc)$type)))
    { type <- as.character(type)[1]
      warning("chart 'type' not specified, assuming \"", type, "\"",
              immediate. = TRUE) }
  if(!exists(paste("stats.", type, sep = ""), mode="function") |
     !exists(paste("sd.", type, sep = ""), mode="function") |
     !exists(paste("limits.", type, sep = ""), mode="function"))
    stop(paste("invalid", type, "control chart. See help(qcc) "))

  if (missing(data.name)) 
     data.name <- deparse(substitute(data))
  data <- data.matrix(data)
  if (missing(sizes)) 
     { if (any(type==c("p", "np", "u")))
          stop(paste("sample 'sizes' must be given for a", type, "Chart"))
       else
          sizes <- apply(data, 1, function(x) sum(!is.na(x)))  }
  else
     { if (length(sizes)==1)
          sizes <- rep(sizes, nrow(data))
       else if (length(sizes) != nrow(data))
                stop("sizes length doesn't match with data") }

  if (missing(labels))
     { if (is.null(rownames(data))) labels <- 1:nrow(data)
       else                         labels <- rownames(data) }

  stats <- paste("stats.", type, sep = "")
  if (!exists(stats, mode="function"))
     stop(paste("function", stats, "is not defined"))
  stats <- do.call(stats, list(data, sizes))
  statistics <- stats$statistics
  if (missing(center)) center <- stats$center

  sd <- paste("sd.", type, sep = "")
  if (!exists(sd, mode="function"))
     stop(paste("function", sd, "is not defined!"))
  missing.std.dev <- missing(std.dev)
  if (missing.std.dev)
     { std.dev <- NULL
       std.dev <- switch(type, 
                         "xbar" = { if(any(sizes > 25)) "RMSDF"
                                    else                "UWAVE-R" },
                         "xbar.one" = "MR",
                         "R" = "UWAVE-R",
                         "S" = "UWAVE-SD",
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

  object <- list(call = call, type = type, 
                 data.name = data.name, data = data, 
                 statistics = statistics, sizes = sizes, 
                 center = center, std.dev = std.dev)
  # check for new data provided and update object
  if (!missing(newdata))
     {   if (missing(newdata.name))
			{newdata.name <- deparse(substitute(newdata))}
       newdata <- data.matrix(newdata)
       if (missing(newsizes))
          { if (any(type==c("p", "np", "u")))
               stop(paste("sample sizes must be given for a", type, "Chart"))
            else
               newsizes <- apply(newdata, 1, function(x) sum(!is.na(x))) }
       else
          { if (length(newsizes)==1)
               newsizes <- rep(newsizes, nrow(newdata))
            else if (length(newsizes) != nrow(newdata))
                     stop("newsizes length doesn't match with newdata") }
       stats <- paste("stats.", type, sep = "")
       if (!exists(stats, mode="function"))
          stop(paste("function", stats, "is not defined"))
       newstats <- do.call(stats, list(newdata, newsizes))$statistics
       if (missing(newlabels))
          { if (is.null(rownames(newdata)))
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

  conf <- nsigmas
  if (!missing(confidence.level))
     conf <- confidence.level
  if (conf >= 1)
     { object$nsigmas <- conf }
  else
     if (conf > 0 & conf < 1)
        { object$confidence.level <- conf } 

  # get control limits
  if (missing(limits))
     { limits <- paste("limits.", type, sep = "")
       if (!exists(limits, mode="function"))
          stop(paste("function", limits, "is not defined"))
       limits <- do.call(limits, list(center = center, std.dev = std.dev,
                                      sizes = sizes, conf = conf)) 
     }
  else 
     { if (!missing.std.dev)
          warning("'std.dev' is not used when limits is given")
       if (!is.numeric(limits))
          stop("'limits' must be a vector of length 2 or a 2-columns matrix")
       limits <- matrix(limits, ncol = 2)
       dimnames(limits) <- list(rep("",nrow(limits)), c("LCL ", "UCL"))
     }

  lcl <- limits[,1]
  ucl <- limits[,2]
  object$limits <- limits
  if (is.function(rules)) violations <- rules(object)
  else                    violations <- NULL
  object$violations <- violations

  class(object) <- "qcc"
  if(plot) plot(object, ...) 
  return(object)
}

print.qcc <- function(x, ...) str(x,1)

summary.qcc <- function(object, digits =  getOption("digits"), ...)
{
  #object <- x   # Argh.  Really want to use 'object' anyway
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
  if(length(center) == 1)
    { cat("\nCenter of group statistics: ", format(center, digits = digits)) }
  else
    { out <- paste(format(center, digits = digits))
      out <- out[which(cumsum(nchar(out)+1) < getOption("width")-40)]      
      out <- paste0(paste(out, collapse = " "), " ...")
      cat("\nCenter of group statistics: ", out, sep = "")
    }
  
  sd <- object$std.dev
  if(length(sd) == 1)
    { cat("\nStandard deviation: ", format(sd, digits = digits), "\n") }
  else
    { out <- paste(format(sd, digits = digits))
      out <- out[which(cumsum(nchar(out)+1) < getOption("width")-40)]
      out <- paste0(paste(out, collapse = " "), " ...")
      cat("\nStandard deviation: ", out, "\n", sep = "")
    }

  newdata.name <- object$newdata.name
  newstats <- object$newstats
  if (!is.null(newstats)) 
     { cat(paste("\nSummary of group statistics in ", 
                 newdata.name, ":\n", sep = ""))
       print(summary(newstats), digits = digits, ...)
       newsizes <- object$newsizes
       if (length(unique(newsizes)) == 1)
          newsizes <- newsizes[1]
       if (length(newsizes) == 1)
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

  limits <- object$limits
  if (!is.null(limits)) 
     { cat("\nControl limits:\n")
       .printShortMatrix(limits, digits = digits, ...) }

  invisible()
}


plot.qcc <- function(x, add.stats = TRUE, chart.all = TRUE, 
                     label.limits = c("LCL ", "UCL"),
                     title, xlab, ylab, ylim, axes.las = 0,
                     digits =  getOption("digits"),
                     restore.par = TRUE, ...) 
{
  object <- x  # Argh.  Really want to use 'object' anyway
  if ((missing(object)) | (!inherits(object, "qcc")))
     stop("an object of class `qcc' is required")

  # collect info from object
  type <- object$type
  std.dev <- object$std.dev
  data.name <- object$data.name
  center <- object$center
  stats <- object$statistics
  limits <- object$limits 
  lcl <- limits[,1]
  ucl <- limits[,2]
  newstats <- object$newstats
  newdata.name <- object$newdata.name
  violations <- object$violations
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
   
  if (missing(title))
     { if (is.null(newstats))
            main.title <- paste(type, "Chart\nfor", data.name)
       else if (chart.all)
                 main.title <- paste(type, "Chart\nfor", data.name, 
                                     "and", newdata.name)
            else main.title <- paste(type, "Chart\nfor", newdata.name) 
     }
  else main.title <- paste(title)

  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  mar <- pmax(oldpar$mar, c(4.1,4.1,3.1,2.1))
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = if(add.stats) pmax(mar, c(7.6,0,0,0)) else mar)

  # plot Shewhart chart
  plot(indices, statistics, type="n",
       ylim = if(!missing(ylim)) ylim 
              else range(statistics, limits, center),
       ylab = if(missing(ylab)) "Group summary statistics" else ylab,
       xlab = if(missing(xlab)) "Group" else xlab, 
       axes = FALSE)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  axis(1, at = indices, las = axes.las,
       labels = if(is.null(names(statistics))) 
                   as.character(indices) else names(statistics))
  axis(2, las = axes.las)
  box()
  top.line <- par("mar")[3]-length(capture.output(cat(main.title)))
  top.line <- top.line - if(chart.all & (!is.null(newstats))) 0.1 else 0.5
  mtext(main.title, side = 3, line = top.line,
        font = par("font.main"), 
        cex  = qcc.options("cex"), 
        col  = par("col.main"))

  lines(indices, statistics, type = "b", pch=20) 

  if(length(center) == 1)
       abline(h = center)
  else lines(indices, center[indices], type="s")

  if(length(lcl) == 1) 
    { abline(h = lcl, lty = 2)
      abline(h = ucl, lty = 2) }
  else 
    { lines(indices, lcl[indices], type="s", lty = 2)
      lines(indices, ucl[indices], type="s", lty = 2) }
  mtext(label.limits, side = 4, at = c(rev(lcl)[1], rev(ucl)[1]), 
        las = 1, line = 0.1, col = gray(0.3), cex = par("cex"))
  mtext("CL", side = 4, at = rev(center)[1], 
        las = 1, line = 0.1, col = gray(0.3), cex = par("cex"))

  if(is.null(qcc.options("violating.runs")))
     stop(".qcc.options$violating.runs undefined. See help(qcc.options).")
  if(length(violations$violating.runs))
    { v <- violations$violating.runs
      if(!chart.all & !is.null(newstats))
        { v <- v - length(stats) 
          v <- v[v>0] }
      points(indices[v], statistics[v], 
             col = qcc.options("violating.runs")$col, 
             pch = qcc.options("violating.runs")$pch) 
    }

  if(is.null(qcc.options("beyond.limits")))
     stop(".qcc.options$beyond.limits undefined. See help(qcc.options).")
  if(length(violations$beyond.limits))
    { v <- violations$beyond.limits
      if(!chart.all & !is.null(newstats))
        { v <- v - length(stats) 
          v <- v[v>0] }
      points(indices[v], statistics[v], 
             col = qcc.options("beyond.limits")$col, 
             pch = qcc.options("beyond.limits")$pch) 
    }

  if(chart.all & (!is.null(newstats)))
    { len.obj.stats <- length(object$statistics)
      len.new.stats <- length(statistics) - len.obj.stats
      abline(v = len.obj.stats + 0.5, lty = 3)
      mtext(# paste("Calibration data in", data.name),
            "Calibration data", cex = par("cex")*0.8,
            at = len.obj.stats/2, line = 0, adj = 0.5)
      mtext(# paste("New data in", object$newdata.name),  
            "New data", cex = par("cex")*0.8, 
            at = len.obj.stats + len.new.stats/2, line = 0, adj = 0.5)
    }

  if(add.stats) 
    { 
      # computes the x margins of the figure region
      plt <- par()$plt; usr <- par()$usr
      px <- diff(usr[1:2])/diff(plt[1:2])
      xfig <- c(usr[1]-px*plt[1], usr[2]+px*(1-plt[2]))
      at.col <- xfig[1] + diff(xfig[1:2])*c(0.10, 0.40, 0.65)
      top.line <- 4.5
      # write info at bottom
      mtext(paste("Number of groups = ", length(statistics), sep = ""), 
            side = 1, line = top.line, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      center <- object$center
      if(length(center) == 1)
        { mtext(paste("Center = ", signif(center[1], digits), sep = ""),
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
      
      if(length(std.dev) == 1)
        { mtext(paste("StdDev = ", signif(std.dev, digits), sep = ""),
                side = 1, line = top.line+2, adj = 0, at = at.col[1],
                font = qcc.options("font.stats"),
                cex = par("cex")*qcc.options("cex.stats"))
        }      
      else
        { mtext("StdDev is variable",
                side = 1, line = top.line+2, adj = 0, at = at.col[1],
                font = qcc.options("font.stats"),
                cex = par("cex")*qcc.options("cex.stats"))
        }
      
      if(length(unique(lcl)) == 1)
        { mtext(paste("LCL = ", signif(lcl[1], digits), sep = ""), 
                side = 1, line = top.line+1, adj = 0, at = at.col[2],
                font = qcc.options("font.stats"),
                cex = par("cex")*qcc.options("cex.stats"))
        }
      else 
        { mtext("LCL is variable", 
                side = 1, line = top.line+1, adj = 0, at = at.col[2],
                font = qcc.options("font.stats"),
                cex = par("cex")*qcc.options("cex.stats"))
        }
      
      if(length(unique(ucl)) == 1)
         { mtext(paste("UCL = ", signif(ucl[1], digits), sep = ""),
                 side = 1, line = top.line+2, adj = 0, at = at.col[2],
                 font = qcc.options("font.stats"),
                 cex = par("cex")*qcc.options("cex.stats")) 
         }
       else 
         { mtext("UCL is variable", 
                 side = 1, line = top.line+2, adj = 0, at = at.col[2],
                 font = qcc.options("font.stats"),
                 cex = par("cex")*qcc.options("cex.stats"))
         }
      
       if(!is.null(violations))
         { mtext(paste("Number beyond limits =",
                       length(unique(violations$beyond.limits))), 
                 side = 1, line = top.line+1, adj = 0, at = at.col[3],
                 font = qcc.options("font.stats"),
                 cex = par("cex")*qcc.options("cex.stats"))
           mtext(paste("Number violating runs =",
                       length(unique(violations$violating.runs))), 
                 side = 1, line = top.line+2, adj = 0, at = at.col[3],
                 font = qcc.options("font.stats"),
                 cex = par("cex")*qcc.options("cex.stats"))
         }
     }

  invisible()
}

#
#  Functions used to compute Shewhart charts statistics
#

.qcc.c4 <- function(n)
{ sqrt(2/(n - 1)) * exp(lgamma(n/2) - lgamma((n - 1)/2)) }
        
# xbar

stats.xbar <- function(data, sizes)
{
  if (missing(sizes))
     sizes <- apply(data, 1, function(x) sum(!is.na(x)))
  statistics <- apply(data, 1, mean, na.rm=TRUE)
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

sd.xbar <- function(data, sizes, std.dev = c("UWAVE-R", "UWAVE-SD", "MVLUE-R", "MVLUE-SD", "RMSDF"))
{
  if (!is.numeric(std.dev))
     std.dev <- match.arg(std.dev)
  if (missing(sizes))
     sizes <- apply(data, 1, function(x) sum(!is.na(x)))
  if (any(sizes == 1))
     stop("group sizes must be larger than one")
  c4 <- .qcc.c4
  if (is.numeric(std.dev))
     { sd <- std.dev }
  else
     { switch(std.dev, 
              "UWAVE-R" = {  R <- apply(data, 1, function(x) 
                                                 diff(range(x, na.rm = TRUE)))
                             d2 <- qcc.options("exp.R.unscaled")[sizes]
                             sd <- sum(R/d2)/length(sizes) }, 
              "UWAVE-SD" = { S <- apply(data, 1, sd, na.rm = TRUE)
                             sd <- sum(S/c4(sizes))/length(sizes) },
              "MVLUE-R"  = { R <- apply(data, 1, function(x) 
                                                 diff(range(x, na.rm = TRUE)))
                             d2 <- qcc.options("exp.R.unscaled")[sizes]
                             d3 <- qcc.options("se.R.unscaled")[sizes]
                             w  <- (d2/d3)^2
                             sd <- sum(R/d2*w)/sum(w) }, 
              "MVLUE-SD" = { S <- apply(data, 1, sd, na.rm = TRUE)
                             w  <- c4(sizes)^2/(1-c4(sizes)^2)
                             sd <- sum(S/c4(sizes)*w)/sum(w) },
              "RMSDF" =    { S <- apply(data, 1, sd, na.rm = TRUE)
                             w  <- sizes-1
                             sd <- sqrt(sum(S^2*w)/sum(w))/c4(sum(w)+1) }
              )
     }
#  if (missing(std.dev))
#     var.within <- apply(data, 1, var, na.rm=TRUE)
#  else 
#     var.within <- std.dev^2
#  var.df <- sum(sizes - 1)
#  if (equal.sd) 
#     { std.dev <- sqrt(sum((sizes - 1) * var.within)/var.df) / c4(var.df + 1) }
#  else 
#     { c <- c4(sizes)/(1 - c4(sizes)^2)
#       std.dev <- sum(c * sqrt(var.within))/sum(c * c4(sizes)) }
  return(sd)
}

limits.xbar <- function(center, std.dev, sizes, conf)
{
  if (length(unique(sizes))==1)
     sizes <- sizes[1]
  se.stats <- std.dev/sqrt(sizes)
  if (conf >= 1) 
     { lcl <- center - conf * se.stats
       ucl <- center + conf * se.stats
     }
  else 
     { if (conf > 0 & conf < 1) 
          { nsigmas <- qnorm(1 - (1 - conf)/2)
            lcl <- center - nsigmas * se.stats
            ucl <- center + nsigmas * se.stats
          }
       else stop("invalid 'conf' argument. See help.")
     }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}

# S chart

stats.S <- function(data, sizes)
{
  if (missing(sizes))
     sizes <- apply(data, 1, function(x) sum(!is.na(x)))
  if(ncol(data)==1) 
    { statistics <- as.vector(data) }
  else 
    { statistics <- sqrt(apply(data, 1, var, na.rm=TRUE)) }
  if (length(sizes == 1))
     sizes <- rep(sizes, length(statistics))
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

sd.S <- function(data, sizes, std.dev = c("UWAVE-SD", "MVLUE-SD", "RMSDF"))
{
  if (!is.numeric(std.dev))
     std.dev <- match.arg(std.dev)
  sd.xbar(data, sizes, std.dev)
}

limits.S <- function(center, std.dev, sizes, conf)
{
  if (length(unique(sizes))==1)
     sizes <- sizes[1]
  c4 <- .qcc.c4
  se.stats <- std.dev * sqrt(1 - c4(sizes)^2)
  if (conf >= 1) 
     { lcl <- pmax(0, center - conf * se.stats)
       ucl <- center + conf * se.stats
     }
  else 
     { if (conf > 0 & conf < 1) 
          { ucl <- std.dev * sqrt(qchisq(1 - (1 - conf)/2, sizes - 1)/
                                  (sizes - 1))
            lcl <- std.dev * sqrt(qchisq((1 - conf)/2, sizes - 1)/
                                  (sizes - 1))
          }
          else stop("invalid conf argument. See help.")
     }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  limits
}

# R Chart 

stats.R <- function(data, sizes)
{
  if (missing(sizes))
     sizes <- apply(data, 1, function(x) sum(!is.na(x)))
  if(ncol(data)==1) 
    { statistics <- as.vector(data) }
  else 
    { statistics <- apply(data, 1, function(x) diff(range(x, na.rm=TRUE))) }
  if (length(sizes == 1))
     sizes <- rep(sizes, length(statistics))
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

sd.R <- function(data, sizes, std.dev = c("UWAVE-R", "MVLUE-R"))
{
  if (!is.numeric(std.dev))
     std.dev <- match.arg(std.dev)
  sd.xbar(data, sizes, std.dev)
}

limits.R <- function(center, std.dev, sizes, conf)
{
  if (length(unique(sizes))==1)
     sizes <- sizes[1]
  se.R.unscaled <- qcc.options("se.R.unscaled")
  Rtab <- length(se.R.unscaled)
  if (conf >= 1) 
     { if (any(sizes > Rtab))
          stop(paste("group size must be less than", 
                      Rtab + 1, "when giving nsigmas"))
       se.R <- se.R.unscaled[sizes] * std.dev
       lcl <- pmax(0, center - conf * se.R)
       ucl <- center + conf * se.R
     }
  else 
     { if (conf > 0 && conf < 1) 
          { ucl <- qtukey(1 - (1 - conf)/2, sizes, 1e100) * std.dev
            lcl <- qtukey((1 - conf)/2, sizes, 1e100) * std.dev
          }
       else stop("invalid conf argument. See help.")
     }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}

# xbar Chart for one-at-time data

stats.xbar.one <- function(data, sizes)
{
  statistics <- as.vector(data)
  center <- mean(statistics)
  list(statistics = statistics, center = center)
}

sd.xbar.one <- function(data, sizes, std.dev = c("MR", "SD"), k = 2)
{
  data <- as.vector(data)
  n <- length(data)
  if(!is.numeric(std.dev)) 
     std.dev <- match.arg(std.dev)
  c4 <- .qcc.c4
  if(is.numeric(std.dev)) 
    { sd <- std.dev }
  else
    { switch(std.dev, 
             "MR" = { d2 <- qcc.options("exp.R.unscaled")
                      if(is.null(d2))
                         stop(".qcc.options$exp.R.unscaled is null")
                      d <- 0
                      for(j in k:n)
                          d <- d+abs(diff(range(data[c(j:(j-k+1))])))
                      sd <- (d/(n-k+1))/d2[k] },
             "SD" = { sd <- sd(data)/c4(n) },
             sd <- NULL)
    }
  return(sd)
}


limits.xbar.one <- function(center, std.dev, sizes, conf)
{
  se.stats <- std.dev
  if (conf >= 1) 
     { lcl <- center - conf * se.stats
       ucl <- center + conf * se.stats
     }
  else 
     { if (conf > 0 & conf < 1) 
          { nsigmas <- qnorm(1 - (1 - conf)/2)
            lcl <- center - nsigmas * se.stats
            ucl <- center + nsigmas * se.stats
          }
       else stop("invalid conf argument. See help.")
     }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}


# p Chart

stats.p <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  list(statistics = data/sizes, center = pbar)
}

sd.p <- function(data, sizes, ...)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  std.dev <- sqrt(pbar * (1 - pbar))
  return(std.dev)
}

limits.p <- function(center, std.dev, sizes, conf)
{ 
  limits.np(center * sizes, std.dev, sizes, conf) / sizes
}

# np Chart

stats.np <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  center <- sizes * pbar
  if (length(unique(center)) == 1)
     center <- center[1]
  list(statistics = data, center = center)
}

sd.np <- function(data, sizes, ...)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  std.dev <- sqrt(sizes * pbar * (1 - pbar))
  if (length(unique(std.dev)) == 1)
     std.dev <- std.dev[1]
  return(std.dev)
}

limits.np <- function(center, std.dev, sizes, conf)
{ 
  sizes <- as.vector(sizes)
  if (length(unique(sizes)) == 1)
     sizes <- sizes[1]
  pbar <- mean(center / sizes)
  if (conf >= 1)
     { tol <- conf * sqrt(pbar * (1 - pbar) * sizes)
       lcl <- pmax(center - tol, 0)
       ucl <- pmin(center + tol, sizes)
     }
  else
     { if (conf > 0 & conf < 1)
          { lcl <- qbinom((1 - conf)/2, sizes, pbar)
            ucl <- qbinom((1 - conf)/2, sizes, pbar, lower.tail = FALSE)
          }
       else stop("invalid conf argument. See help.")
     }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}

# c Chart

stats.c <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  if (length(unique(sizes)) != 1)
     stop("all sizes must be be equal for a c chart")
  statistics <- data
  center <- mean(statistics)
  list(statistics = statistics, center = center)
}

sd.c <- function(data, sizes, ...)
{
  data <- as.vector(data)
  std.dev <- sqrt(mean(data))
  return(std.dev)
}

limits.c <- function(center, std.dev, sizes, conf)
{
  if (conf >= 1) 
     { lcl <- center - conf * sqrt(center)
       lcl[lcl < 0] <- 0
       ucl <- center + conf * sqrt(center)
     }
  else 
     { if (conf > 0 & conf < 1) 
          { ucl <- qpois(1 - (1 - conf)/2, center)
            lcl <- qpois((1 - conf)/2, center)
          }
       else stop("invalid conf argument. See help.")
     }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}

# u Chart

stats.u <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  statistics <- data/sizes
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

sd.u <- function(data, sizes, ...)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  std.dev <- sqrt(sum(data)/sum(sizes))
  return(std.dev)
}

limits.u <- function(center, std.dev, sizes, conf)
{
  sizes <- as.vector(sizes)
  if (length(unique(sizes))==1)
     sizes <- sizes[1]
  limits.c(center * sizes, std.dev, sizes, conf) / sizes
}

#
# Functions used to signal points out of control 
#

shewhart.rules <- function(object, limits = object$limits, run.length = qcc.options("run.length"))
{
# Return a list of cases beyond limits and violating runs
  bl <- beyond.limits(object, limits = limits)
  vr <- violating.runs(object, run.length = run.length)
  list(beyond.limits = bl, violating.runs = vr)
}

beyond.limits <- function(object, limits = object$limits)
{
# Return cases beyond limits
  statistics <- c(object$statistics, object$newstats) 
  lcl <- limits[,1]
  ucl <- limits[,2]
  index.above.ucl <- seq(along = statistics)[statistics > ucl]
  index.below.lcl <- seq(along = statistics)[statistics < lcl]
  return(c(index.above.ucl, index.below.lcl))
}

violating.runs <- function(object, run.length = qcc.options("run.length"))
{
# Return indices of points violating runs
  if(run.length == 0)
    return(numeric())
  center <- object$center
  statistics <- c(object$statistics, object$newstats)
  cl <- object$limits
  diffs <- statistics - center
  diffs[diffs > 0] <- 1
  diffs[diffs < 0] <- -1
  runs <- rle(diffs)
  vruns <- rep(runs$lengths >= run.length, runs$lengths)
  vruns.above <- (vruns & (diffs > 0))
  vruns.below <- (vruns & (diffs < 0))
  rvruns.above <- rle(vruns.above)
  rvruns.below <- rle(vruns.below)
  vbeg.above <- cumsum(rvruns.above$lengths)[rvruns.above$values] -
                      (rvruns.above$lengths - run.length)[rvruns.above$values]
  vend.above <- cumsum(rvruns.above$lengths)[rvruns.above$values]
  vbeg.below <- cumsum(rvruns.below$lengths)[rvruns.below$values] -
                      (rvruns.below$lengths - run.length)[rvruns.below$values]
  vend.below <- cumsum(rvruns.below$lengths)[rvruns.below$values]
  violators <- numeric()
  if (length(vbeg.above)) 
     { for (i in 1:length(vbeg.above))
           violators <- c(violators, vbeg.above[i]:vend.above[i]) }
  if (length(vbeg.below)) 
     { for (i in 1:length(vbeg.below))
           violators <- c(violators, vbeg.below[i]:vend.below[i]) }
  return(violators)
}

#-------------------------------------------------------------------#
#                                                                   #
#          Operating Characteristic Function                        #
#                                                                   #
#-------------------------------------------------------------------#

oc.curves <- function(object, ...)
{
# Draws the operating characteristic curves for the qcc object 

  if ((missing(object)) | (!inherits(object, "qcc")))
     stop("an object of class 'qcc' is required")

  size <- unique(object$sizes)
  if (length(size)>1)
     stop("Operating characteristic curves available only for equal sample sizes!")

  beta <- switch(object$type,
                 xbar = oc.curves.xbar(object, ...),
                 R    = oc.curves.R(object, ...),
                 S    = oc.curves.S(object, ...),
                 np   =,
                 p    = oc.curves.p(object, ...),
                 u    =,
                 c    = oc.curves.c(object, ...))
  if (is.null(beta))
     stop("Operating characteristic curves not available for this type of chart.")

  invisible(beta)
}

oc.curves.xbar <- function(object, n, c = seq(0, 5, length=101), nsigmas = object$nsigmas, identify=FALSE, restore.par=TRUE)
{
# Draw the operating-characteristic curves for the xbar-chart with nsigmas
# limits. The values on the vertical axis give the probability of not detecting
# a shift of c*sigma in the mean on the first sample following the shift.

  if (!(object$type=="xbar"))
     stop("not a `qcc' object of type \"xbar\".")

  size <- unique(object$sizes)
  if (length(size) > 1)
     stop("Operating characteristic curves available only for equal sample sizes!")
  if (missing(n))
     n <- unique(c(size, c(1,5,10,15,20)))
  if (is.null(nsigmas))
     nsigmas <- qnorm(1 - (1 - object$confidence.level) / 2)

  beta <- matrix(as.double(NA), length(n), length(c))
  for (i in 1:length(n))
      beta[i,] <- pnorm(nsigmas-c*sqrt(n[i])) - pnorm(-nsigmas-c*sqrt(n[i]))
  rownames(beta) <- paste("n=",n,sep="")
  colnames(beta) <- c

  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = pmax(oldpar$mar, c(4.1,4.1,2.1,2.1)))

  plot(c, beta[1,], type="n",
       ylim = c(0,1), xlim = c(0,max(c)),
       xlab = "Process shift (std.dev)",
       ylab = "Prob. type II error ")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  box()
  mtext(paste("OC curves for", object$type, "Chart"), 
        side = 3, line = par("mar")[3]/3,
        font = par("font.main"), 
        cex  = qcc.options("cex"), 
        col  = par("col.main"))
  for(i in 1:length(n))
     lines(c, beta[i,], type = "l", lty=i)
  beta <- t(beta)
  names(dimnames(beta)) <- c("shift (std.dev)", "sample size")

  if (identify)
     { cs <- rep(c,length(n))
       betas <- as.vector(beta)
       labels <- paste("c=", formatC(cs, 2, flag="-"), 
                       ": beta=", formatC(betas, 4, flag="-"), 
                       ", ARL=", formatC(1/(1-betas), 2, flag="-"), sep="")
       i <- identify(cs, betas, labels, pos=4, offset=0.2)
       apply(as.matrix(labels[i$ind]), 1, cat, "\n")
     }
  else
     { legend(max(c), 1, legend = paste("n =", n), 
              bg = qcc.options("bg.figure"),
              lty = 1:length(n), xjust = 1, yjust = 1)
     }
  invisible(beta)
}

oc.curves.p <- function(object, nsigmas = object$nsigmas, identify=FALSE, restore.par=TRUE)
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
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = pmax(oldpar$mar, c(4.1,4.1,2.1,2.1)))

  plot(p, beta, type = "n", 
       ylim = c(0,1), xlim = c(0,1),
       xlab = expression(p), 
       ylab = "Prob. type II error ")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  box()
  mtext(paste("OC curves for", object$type, "Chart"), 
        side = 3, line = par("mar")[3]/3,
        font = par("font.main"), 
        cex  = qcc.options("cex"), 
        col  = par("col.main"))

  lines(p, beta)
  lines(rep(p[which.max(beta)], 2), c(0, max(beta)), lty = 2)

  warning("Some computed values for the type II error have been rounded due to the discreteness of the binomial distribution. Thus, some ARL values might be meaningless.")

  if (identify)
     { labels <- paste("p=", formatC(p, 2, flag="-"), 
                       ": beta=", formatC(beta, 4, flag="-"), 
                       ", ARL=", formatC(1/(1-beta), 2, flag="-"), sep="")
       i <- identify(p, beta, labels, pos=4, offset=0.2)
       apply(as.matrix(labels[i$ind]), 1, cat, "\n")
     }
  invisible(beta)  
}

oc.curves.c <- function(object, nsigmas = object$nsigmas, identify=FALSE, restore.par=TRUE)
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
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = pmax(oldpar$mar, c(4.1,4.1,2.1,2.1)))

  plot(lambda, beta, type = "n", 
       ylim = c(0,1), xlim = range(lambda),
       xlab = "Mean", 
       ylab = "Prob. type II error ")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  box()
  mtext(paste("OC curves for", object$type, "Chart"), 
        side = 3, line = par("mar")[3]/3,
        font = par("font.main"), 
        cex  = qcc.options("cex"), 
        col  = par("col.main"))

  lines(lambda, beta)     
  lines(rep(lambda[which.max(beta)], 2), c(0, max(beta)), lty = 2)

  warning("Some computed values for the type II error have been rounded due to the discreteness of the Poisson distribution. Thus, some ARL values might be meaningless.")

  if (identify)
     { labels <- paste("lambda=", formatC(lambda, 0, flag="-"), 
                       ": beta=", formatC(beta, 4, flag="-"), 
                       ", ARL=", formatC(1/(1-beta), 2, flag="-"), sep="")
       i <- identify(lambda, beta, labels, pos=4, offset=0.2)
       apply(as.matrix(labels[i$ind]), 1, cat, "\n")
     }
  invisible(beta)
}

oc.curves.R <-
function(object, n, c = seq(1, 6, length=101), nsigmas = object$nsigmas, identify=FALSE, restore.par=TRUE)
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
     n <- unique(c(size, c(2,5,10,15,20)))
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

  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = pmax(oldpar$mar, c(4.1,4.1,2.1,2.1)))

  plot(c, beta[,1], type="n",
       ylim = c(0,1), xlim = c(1,max(c)),
       xlab = "Process scale multiplier",
       ylab = "Prob. type II error ")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
       col = qcc.options("bg.figure"))
  box()
  mtext(paste("OC curves for", object$type, "Chart"), 
        side = 3, line = par("mar")[3]/3,
        font = par("font.main"), 
        cex  = qcc.options("cex"), 
        col  = par("col.main"))
  
  matlines(c, beta, lty = 1:length(n), col = 1)

  names(dimnames(beta)) <- c("scale multiplier", "sample size")

  if (identify)
     { cs <- rep(c,length(n))
       betas <- as.vector(beta)
       labels <- paste("c=", formatC(cs, 2, flag="-"),
                       ": beta=", formatC(betas, 4, flag="-"),
                       ", ARL=", formatC(1/(1-betas), 2, flag="-"), sep="")
       i <- identify(cs, betas, labels, pos=4, offset=0.2)
       apply(as.matrix(labels[i$ind]), 1, cat, "\n")
     }
  else
     { legend(max(c), 1, legend = paste("n =", n),
              bg = qcc.options("bg.figure"),
              lty = 1:length(n), xjust = 1, yjust = 1)
     }
  invisible(beta)
}

oc.curves.S <- function(object, n, c = seq(1, 6, length=101), nsigmas = object$nsigmas, identify=FALSE, restore.par=TRUE)
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
     n <- unique(c(size, c(2,5,10,15,20)))
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
    { c4 <- .qcc.c4
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

  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = pmax(oldpar$mar, c(4.1,4.1,2.1,2.1)))

  plot(c, beta[,1], type="n",
       ylim = c(0,1), xlim = c(1,max(c)),
       xlab = "Process scale multiplier",
       ylab = "Prob. type II error ")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
       col = qcc.options("bg.figure"))
  box()
  mtext(paste("OC curves for", object$type, "Chart"), 
        side = 3, line = par("mar")[3]/3,
        font = par("font.main"), 
        cex  = qcc.options("cex"), 
        col  = par("col.main"))
  matlines(c, beta, lty = 1:length(n), col = 1)

  names(dimnames(beta)) <- c("scale multiplier", "sample size")

  if (identify)
     { cs <- rep(c,length(n))
       betas <- as.vector(beta)
       labels <- paste("c=", formatC(cs, 2, flag="-"),
                       ": beta=", formatC(betas, 4, flag="-"),
                       ", ARL=", formatC(1/(1-betas), 2, flag="-"), sep="")
       i <- identify(cs, betas, labels, pos=4, offset=0.2)
       apply(as.matrix(labels[i$ind]), 1, cat, "\n")
     }
  else
     { legend(max(c), 1, legend = paste("n =", n),
              bg = qcc.options("bg.figure"),
              lty = 1:length(n), xjust = 1, yjust = 1)
     }
  invisible(beta)
}

#-------------------------------------------------------------------#
#                                                                   #
# Miscellaneous functions                                           #
#                                                                   #
#-------------------------------------------------------------------#

qcc.groups <- function(data, sample)
{
  if(length(data)!=length(sample))
    stop("data and sample must be vectors of equal length")
  x <- lapply(split(data, sample), as.vector)
  lx <- sapply(x, length)
  for(i in which(lx != max(lx)))
      x[[i]] <- c(x[[i]], rep(NA, max(lx)-lx[i]))
  x <- t(sapply(x, as.vector))
  return(x)
}

qcc.overdispersion.test <- function(x, size, 
                            type=ifelse(missing(size), "poisson", "binomial"))
{
  type <- match.arg(type, c("poisson", "binomial"))
  if (type=="binomial" & missing(size))
     stop("binomial data require argument \"size\"")
  if (!missing(size))
     if (length(x) != length(size))   
        stop("arguments \"x\" and \"size\" must be vector of same length")

  n <- length(x)
  obs.var <- var(x)
  if (type=="binomial")
     { p <- sum(x)/sum(size)
       theor.var <- mean(size)*p*(1-p) }
  else if (type=="poisson")
          { theor.var <- mean(x) }
       else
          stop("invalid \"type\" argument. See help.")

  D <- (obs.var * (n-1)) / theor.var
  p.value <- 1-pchisq(D, n-1)

  out <- matrix(c(obs.var/theor.var, D, signif(p.value,5)), 1, 3)
  rownames(out) <- paste(type, "data")
  colnames(out) <- c("Obs.Var/Theor.Var", "Statistic", "p-value") 
  names(dimnames(out)) <- c(paste("Overdispersion test"), "")
  return(out)
}


