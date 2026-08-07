#' Quality Control Charts
#'
#' Create an object of class `'qcc'` to perform statistical quality
#' control. This object may then be used to plot Shewhart charts, drawing OC
#' curves, computes capability indices, and more.
#'
#' Numeric `rules` values are interpreted within `rule.set`. By
#' default, `rules = c(1,4)` applies Western Electric rules 1 and 4 for
#' backward compatibility. Nelson rules can be requested with `rules =
#' 1:8, rule.set = "nelson"`.
#'
#' @export
#' @inheritParams chart_common data newdata newsizes center
#' @param type a character string specifying the group statistics to compute.
#' Available methods are:
#'
#' | Type | Statistic charted | Chart description |
#' | --- | --- | --- |
#' | `"xbar"` | mean | means of a continuous process variable |
#' | `"R"` | range | ranges of a continuous process variable |
#' | `"S"` | standard deviation | standard deviations of a continuous variable |
#' | `"xbar.one"` | mean | one-at-time data of a continuous process variable |
#' | `"p"` | proportion | proportion of nonconforming units |
#' | `"np"` | count | number of nonconforming units |
#' | `"c"` | count | nonconformities per unit |
#' | `"u"` | count | average nonconformities per unit |
#' | `"g"` | count | number of non-events between events |
#'
#' Furthermore, a user specified type of chart, say `"newchart"`, can be
#' provided. This requires the definition of `"stats.newchart"`,
#' `"sd.newchart"`, and `"limits.newchart"`. As an example, see
#' [stats.xbar()].
#' @param sizes a value or a vector of values specifying the sample sizes
#' associated with each group. For continuous data provided as data frame or
#' matrix the sample sizes are obtained counting the non-`NA` elements of
#' each row. For `"p"`, `"np"` and `"u"` charts the argument
#' `sizes` is required.
#' @param std.dev a value or an available method specifying the within-group
#' standard deviation(s) of the process. Several methods are available for
#' estimating the standard deviation in case of a continuous process variable;
#' see [sd.xbar()], [sd.xbar.one()], [sd.R()], and [sd.S()].
#' @param limits a two-values vector specifying control limits.
#' @param nsigmas a numeric value specifying the number of sigmas to use for
#' computing control limits. It is ignored when the `confidence.level`
#' argument is provided.
#' @param confidence.level a numeric value between 0 and 1 specifying the
#' confidence level of the computed probability limits.
#' @param rules a value or a vector of values specifying the rules to apply to
#' the chart. See [qccRules()] for possible values and their meaning.
#' @param rule.set a character string specifying how numeric `rules`
#' values are interpreted. The default is `"western-electric"` specifying
#' Western Electric rules 1 through 4. Use `"nelson"` to apply Nelson
#' rules 1 through 8.
#' @param label.center a character specifying the label for center line.
#' @param label.limits a character vector specifying the labels for control
#' limits.
#' @param x an object of class `'qcc'`.
#' @param object an object of class `'qcc'`.
#' @param ... additional arguments to be passed to the generic function.
#' @return Returns an object of class `'qcc'`.
#' @author Luca Scrucca
#' @family control charts
#' @seealso [qccRules()], [ocCurves()], [processCapability()], [qccGroups()]
#' @references `r refs("mason_young_2002", "montgomery2013", "ryan_2011", "scrucca_2004", "wetherill_brown_1991")`
#' @examples
#' ##  Continuous data 
#' diameter  = qccGroups(data = pistonrings, diameter, sample)
#'
#' (q  = qcc(diameter[1:25,], type="xbar"))
#' plot(q)
#'
#' (q  = qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,]))
#' plot(q)
#'
#' q = qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,])
#' plot(q, chart.all=FALSE)
#'
#' plot(qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,], nsigmas=2))
#'
#' plot(qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,], confidence.level=0.99))
#'
#' q <- qcc(diameter[1:25,], type="R")
#' q
#' plot(q)
#'
#' plot(qcc(diameter[1:25,], type="R", newdata=diameter[26:40,]))
#'
#' plot(qcc(diameter[1:25,], type="S"))
#'
#' plot(qcc(diameter[1:25,], type="S", newdata=diameter[26:40,]))
#'
#' plot(qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,], rules = 1:4))
#'
#' # variable control limits
#' out  = c(9, 10, 30, 35, 45, 64, 65, 74, 75, 85, 99, 100)
#' diameter  = qccGroups(data = pistonrings[-out,], diameter, sample)
#' plot(qcc(diameter[1:25,], type="xbar"))
#' plot(qcc(diameter[1:25,], type="R"))
#' plot(qcc(diameter[1:25,], type="S"))
#' plot(qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,]))
#' plot(qcc(diameter[1:25,], type="R", newdata=diameter[26:40,]))
#' plot(qcc(diameter[1:25,], type="S", newdata=diameter[26:40,]))
#'
#' # to customize a Shewhart chart use
#' q <- qcc(diameter[1:25, ], type = "xbar")
#' graph = plot(q, ylim = c(73.9, 74.1)) 
#' # returned object is of class "patchwork", then add geom_* layer to the first element
#' graph[[1]] <- graph[[1]] + 
#'   geom_hline(yintercept = c(73.95, 74.05), lty = 2) 
#' graph
#'
#' ##  Attribute data 
#' q <- with(orangejuice, qcc(D[trial], sizes=size[trial], type="p"))
#' q
#' plot(q)
#'
#' # remove out-of-control points (see help(orangejuice) for the reasons)
#' outofctrl <- c(15,23)
#' q1 <- with(orangejuice[-outofctrl,], 
#'            qcc(D[trial], sizes=size[trial], type="p"))
#' plot(q1)
#' q1 <- with(orangejuice[-outofctrl,], 
#'            qcc(D[trial], sizes=size[trial], type="p",
#'                newdata=D[!trial], newsizes=size[!trial]))
#' plot(q1)
#'
#' data(orangejuice2)
#' q2  <- with(orangejuice2, 
#'            qcc(D[trial], sizes=size[trial], type="p"))
#' plot(q2)
#' q2  <- with(orangejuice2, 
#'            qcc(D[trial], sizes=size[trial], type="p", 
#'                newdata=D[!trial], newsizes=size[!trial]))
#' plot(q2)
#'
#' data(circuit)
#' plot(with(circuit, qcc(x[trial], sizes=size[trial], type="c")))
#'
#' # remove out-of-control points (see help(circuit) for the reasons)
#' outofctrl  = c(15,23)
#' q1  <- with(orangejuice[-outofctrl,], 
#'            qcc(D[trial], sizes=size[trial], type="p"))
#' plot(q1)
#' q1  <- with(orangejuice[-outofctrl,], 
#'            qcc(D[trial], sizes=size[trial], type="p",
#'                newdata=D[!trial], newsizes=size[!trial]))
#' plot(q1)
#'
#' outofctrl  = c(6,20)
#' q1  <- with(circuit[-outofctrl,], 
#'            qcc(x[trial], sizes=size[trial], type="c"))
#' plot(q1)
#' q1  <- with(circuit[-outofctrl,], 
#'            qcc(x[trial], sizes=size[trial], type="c", 
#'                newdata = x[!trial], newsizes = size[!trial]))
#' plot(q1)
#' q1  <- with(circuit[-outofctrl,], 
#'            qcc(x[trial], sizes=size[trial], type="u", 
#'            newdata = x[!trial], newsizes = size[!trial]))
#' plot(q1)
#'
#' data(pcmanufact)
#' q1  <- with(pcmanufact, qcc(x, sizes=size, type="u"))
#' q1
#' plot(q1)
#'
#' data(dyedcloth)
#' # variable control limits
#' plot(with(dyedcloth, qcc(x, sizes=size, type="u")))
#' # standardized control chart
#' q  <- with(dyedcloth, qcc(x, sizes=size, type="u"))
#' z  <- (q$statistics - q$center)/sqrt(q$center/q$size)
#' plot(qcc(z, sizes = 1, type = "u", center = 0, std.dev = 1, limits = c(-3,3)),
#'      title = "Standardized u chart")
#'     
#' ##  Continuous one-at-time data 
#'
#' q <- with(viscosity, qcc(viscosity[trial], type = "xbar.one"))
#' q
#' plot(q)
#' # batch 4 is out-of-control because of a process temperature controller
#' # failure; remove it and recompute
#' viscosity  <- viscosity[-4,]
#' plot(with(viscosity, 
#'           qcc(viscosity[trial], type = "xbar.one", newdata = viscosity[!trial])))
qcc <- function(data, 
                type = c("xbar", "R", "S", "xbar.one", 
                         "p", "np", "c", "u", "g"),
                sizes, center, std.dev, limits, 
                newdata, newsizes, 
                nsigmas = 3, confidence.level, 
                rules = c(1,4),
                rule.set = c("western-electric", "nelson"), ...)
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

  data.name <- deparse(substitute(data))
  data <- data.matrix(data)
  if (missing(sizes)) 
     { if (any(type==c("p", "np", "u")))
          stop(paste("sample 'sizes' must be given for a", type, "Chart"))
       else
          sizes <- as.integer(rowSums(!is.na(data))) }
  else
     { if (length(sizes)==1)
          sizes <- rep(sizes, nrow(data))
       else if (length(sizes) != nrow(data))
                stop("sizes length doesn't match with data") }

  labels <- rownames(data) %||% 1:nrow(data)

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

  stopifnot(length(labels) == length(statistics))
  names(statistics) <-  rownames(data) <-  labels
  names(dimnames(data)) <- list("Group", "Samples")
  rules <- as.numeric(rules)
  rule.set <- match.arg(rule.set)
  
  # create object of class 'qcc'
  object <- list(call = call, type = type,
                 data.name = data.name, data = data, 
                 statistics = statistics, sizes = sizes, 
                 center = center, std.dev = std.dev,
                 rules = rules, rule.set = rule.set)
  class(object) <- "qcc"

  # check for new data provided and update object
  if(!missing(newdata))
  { 
    newdata.name <- deparse(substitute(newdata))
    newdata <- data.matrix(newdata)
    if(missing(newsizes))
    { 
      if(any(type==c("p", "np", "u")))
        stop(paste("sample 'newsizes' must be given for a", type, "Chart"))
      else
        newsizes <- as.integer(rowSums(!is.na(newdata)))
    } else
    { 
      if(length(newsizes)==1)
        newsizes <- rep(newsizes, nrow(newdata))
      else 
        if(length(newsizes) != nrow(newdata))
         stop("newsizes length doesn't match with newdata") 
    }
    stats <- paste("stats.", type, sep = "")
    if(!exists(stats, mode="function"))
      stop(paste("function", stats, "is not defined"))
    newstats <- do.call(stats, list(newdata, newsizes))$statistics
    if(is.null(rownames(newdata)))
    { 
      start <- length(statistics)
      newlabels <- seq(start+1, start+length(newstats)) 
    } else
    { 
      newlabels <- rownames(newdata)
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
  
  conf <- nsigmas
  if (!missing(confidence.level))
     conf <- confidence.level
  if (conf >= 1)
     { object$nsigmas <- conf }
  else
     if (conf > 0 & conf < 1)
        { object$confidence.level <- conf
          object$nsigmas <- qnorm(1 - (1 - conf)/2) }

  # get control limits
  if (missing(limits))
     { limits <- paste("limits.", type, sep = "")
       if (!exists(limits, mode="function"))
          stop(paste("function", limits, "is not defined"))
       limits <- do.call(limits, list(center = center, 
                                      std.dev = std.dev,
                                      sizes = sizes, 
                                      nsigmas = object$nsigmas,
                                      conf = object$confidence.level)) 
     }
  else 
     { if (!missing.std.dev)
          warning("'std.dev' is not used when limits is given")
       if (!is.numeric(limits))
          stop("'limits' must be a vector of length 2 or a 2-columns matrix")
       limits <- new_limits(matrix(limits, ncol = 2))
     }
  object$limits <- limits
  
  # identify violating rules observations
  object$violations <- qccRules(object)

  return(object)
}


#' @rdname qcc
#' @export
#' @export print.qcc
print.qcc <- function(x, digits = getOption("digits"), ...)
{
  object <- x   # Argh.  Really want to use 'object' anyway
  cat(cli::rule(left = crayon::bold("Quality Control Chart"), 
                width = min(getOption("width"),50)), "\n\n")

  data.name <- object$data.name
  type <- object$type
  statistics <- object$statistics
  cat("Chart type                 =", type, "\n")
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
  if(length(center) == 1)
  { 
    cat("Center of group statistics =", signif(center, digits = digits), "\n")   } else
  { 
    ng  <- length(center)
    center <- paste(signif(center, digits = digits))
    out <- if(ng > 3) c(center[1:2], "...", center[ng]) else center
    out <- out[which(cumsum(nchar(out)+1) < getOption("width")-40)] 
    out <-paste(out, collapse = " ")
    cat("Center of group statistics = ", out, "\n", sep = "")
  }
  
  sd <- object$std.dev
  if(length(sd) == 1)
  { 
    cat("Standard deviation         =", signif(sd, digits = digits), "\n") 
  } else
  { 
    ng  <- length(sd)
    sd <- paste(signif(sd, digits = digits))
    out <- if(ng > 3) c(sd[1:2], "...", sd[ng]) else sd
    out <- out[which(cumsum(nchar(out)+1) < getOption("width")-40)] 
    out <-paste(out, collapse = " ")
    cat("Standard deviation         = ", out, "\n", sep = "")
  }

  newdata.name <- object$newdata.name
  newstats <- object$newstats
  if (!is.null(newstats)) 
  { 
    cat("\nNew data (phase II)        =", newdata.name, "\n")
    cat("Number of groups           =", length(newstats), "\n")
    newsizes <- object$newsizes
    if (length(unique(newsizes)) == 1)
      newsizes <- newsizes[1]
    if (length(newsizes) == 1)
    {
      cat("Group sample size          =", signif(newsizes), "\n")
    } else 
    { cat("Group sample sizes         =")
      new.tab <- table(newsizes)
      print(matrix(c(as.numeric(names(new.tab)), new.tab),
                   ncol = length(new.tab), byrow = TRUE, 
                   dimnames = list(c("  sizes", "  counts"),
                                   character(length(new.tab)))), 
            digits = digits, ...)
    }
  }
  
  limits <- object$limits
  if(!is.null(limits)) 
  { 
    cat("\nControl limits at nsigmas  =", object$nsigmas, "\n")    
    .printShortMatrix(limits, digits = digits, ...) 
  }

  invisible()
}

#' @rdname qcc
#' @export
#' @export summary.qcc
summary.qcc <- function(object, ...) print.qcc(object, ...)


#' @rdname qcc
#' @export
#' @export plot.qcc
#' @inheritParams plot_common 
plot.qcc <- function(x, xtime = NULL,
                     add.stats = qcc.options("add.stats"), 
                     chart.all = qcc.options("chart.all"), 
                     fill = qcc.options("fill"),
                     label.center = "CL",
                     label.limits = c("LCL ", "UCL"), 
                     title, xlab, ylab, xlim, ylim,
                     digits = getOption("digits"),
                     ...) 
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
  rules <- object$rules
  rule.set <- object$rule.set %||% "western-electric"
  rule.set <- match.arg(rule.set, c("western-electric", "nelson"))
  statistics <- c(stats, newstats)
  plot.index <- qcc_plot_index(
    n_phase1 = length(stats),
    n_phase2 = length(newstats),
    xtime = xtime,
    chart_all = chart.all
  )
  groups <- plot.index$group
  
  if(missing(title))
  { 
    if(is.null(newstats))
      title <- paste(type, "chart for", data.name)
    else if(chart.all)
           title <- paste(type, "chart for", data.name, "and", newdata.name)
         else 
           title <- paste(type, "chart for", newdata.name) 
  }
  
  violation.values <- ifelse(is.na(violations), 0, violations)
  violation.levels <- sort(unique(c(0, violation.values)))
  rule.options <- qcc.options("rules")
  colour.values <- setNames(
    c("black", rule.options$col),
    c("0", seq_along(rule.options$col))
  )
  shape.values <- setNames(
    c(20, rule.options$pch),
    c("0", seq_along(rule.options$pch))
  )

  plot.data <- data.frame(
    stat = statistics,
    center = center,
    lcl = lcl,
    ucl = ucl,
    violations = factor(violation.values, levels = violation.levels)
  )
  df <- cbind(
    plot.index,
    plot.data[plot.index$row, , drop = FALSE]
  )
  
  if(missing(ylim))
    ylim <- extendrange(c(df$stat, df$lcl, df$ucl))
  if(missing(xlim))
    xlim <- extendrange(df$group)

  plot <- 
    ggplot(data = df, aes(x = .data[["group"]], 
                          y = .data[["stat"]])) +
    geom_line() +
    geom_point(aes(colour = .data[["violations"]], 
                   shape = .data[["violations"]]), 
               size = 2) +
    scale_colour_manual(values = colour.values,
                        breaks = levels(df$violations)) +
    scale_shape_manual(values = shape.values,
                       breaks = levels(df$violations)) +
    labs(title = title, subtitle = "",
         x = if(missing(xlab)) "Group" else xlab,
         y = if(missing(ylab)) "Group summary statistics" else ylab) +
    coord_cartesian(xlim = xlim, ylim = ylim,
                    expand = FALSE, clip = "off") +
    theme_qcc()

  plot <- plot + scale_x_qcc(groups, xlim)
        
  # draw control limits
  has.rule <- function(x) any(rules %in% x)

  if(has.rule(1))
  { 
    dx <- min(diff(df$group))/2
    x1 <- x2 <- c(xlim[1], df$group[-length(df$group)]+dx, xlim[2])
    y1 <- c(df$lcl, df$lcl[length(df$group)])
    y2 <- c(df$ucl, df$ucl[length(df$group)])
    xp1 <- rep(x1, each=2)[-1]
    xp2 <- rep(x2, each=2)[-1]
    yp1 <- rep(y1, each=2)[-2*length(y1)]
    yp2 <- rep(y2, each=2)[-2*length(y2)]
    if(fill)
    { 
      # fill the in-control area
      plot <- plot + 
        geom_polygon(data = data.frame(x = c(xp1,rev(xp2)), 
                                       y = c(yp1,rev(yp2))),
                     aes(x = .data[["x"]], 
                         y = .data[["y"]]), 
                     fill = adjustcolor(qcc.options("zones")$fill, alpha.f=0.2),
                     col = NA)
    } else
    {
      plot <- plot + 
        geom_step(data = data.frame(x = x1, 
                                    y = y1),
                  aes(x = .data[["x"]], 
                      y = .data[["y"]]), 
                  lty = qcc.options("zones")$lty[1],
                  col = qcc.options("zones")$col[1])
      plot <- plot + 
        geom_step(data = data.frame(x = x2, 
                                    y = y2),
                  aes(x = .data[["x"]], 
                      y = .data[["y"]]), 
                  lty = qcc.options("zones")$lty[1],
                  col = qcc.options("zones")$col[1])
    }

    plot <- plot + 
      annotate("text", x = Inf, 
               y = c(rev(df$center)[1], rev(df$lcl)[1], rev(df$ucl)[1]),
               label = c(label.center, label.limits),
               col = gray(0.3), size = 10 * 5/14,
               hjust = -0.2, vjust = 0.5)
  }
  
  # draw 2-sigma warning limits
  if(if(rule.set == "nelson") has.rule(5) else has.rule(2))
  { 
    limits.2sigma <- do.call(paste("limits.", object$type, sep = ""), 
                             list(center = object$center, 
                                  std.dev = object$std.dev,
                                  sizes = c(object$sizes, object$newsizes), 
                                  nsigmas = object$nsigmas*2/3))
    dx <- min(diff(df$group))/2
    x1 <- x2 <- c(xlim[1], df$group[-length(df$group)]+dx, xlim[2])
    if(nrow(limits.2sigma)==1)
    {
      y1 <- rep(limits.2sigma[1,1], length(df$group)+1)
      y2 <- rep(limits.2sigma[1,2], length(df$group)+1)
    } else
    {
      limits.2sigma <- limits.2sigma[plot.index$row, , drop = FALSE]
      y1 <- c(limits.2sigma[,1], limits.2sigma[length(df$group),1])
      y2 <- c(limits.2sigma[,2], limits.2sigma[length(df$group),2])
    }
    xp1 <- rep(x1, each=2)[-1]
    xp2 <- rep(x2, each=2)[-1]
    yp1 <- rep(y1, each=2)[-2*length(y1)]
    yp2 <- rep(y2, each=2)[-2*length(y2)]
    if(fill)
    { # fill the in-control area
      plot <- plot + 
        geom_polygon(data = data.frame(x = c(xp1,rev(xp2)), 
                                       y = c(yp1,rev(yp2))),
                     aes(x = .data[["x"]], 
                         y = .data[["y"]]), 
                     fill = adjustcolor(qcc.options("zones")$fill, alpha.f=0.2),
                     col = NA)
    } else
    {
      plot <- plot + 
        geom_step(data = data.frame(x = x1, y = y1),
                           aes(x = .data[["x"]], 
                               y = .data[["y"]]), 
                           lty = qcc.options("zones")$lty[2],
                           col = qcc.options("zones")$col[2])
      plot <- plot + 
        geom_step(data = data.frame(x = x2, y = y2),
                           aes(x = .data[["x"]], 
                               y = .data[["y"]]), 
                           lty = qcc.options("zones")$lty[2],
                           col = qcc.options("zones")$col[2])
    }
  }
  
  # draw 1-sigma warning limits
  if(if(rule.set == "nelson") has.rule(c(6, 7, 8)) else has.rule(3))
  { 
    limits.2sigma <- do.call(paste("limits.", object$type, sep = ""), 
                             list(center = object$center, 
                                  std.dev = object$std.dev,
                                  sizes = c(object$sizes, object$newsizes), 
                                  nsigmas = object$nsigmas*1/3))
    dx <- min(diff(df$group))/2
    x1 <- x2 <- c(xlim[1], df$group[-length(df$group)]+dx, xlim[2])
    if(nrow(limits.2sigma)==1)
    {
      y1 <- rep(limits.2sigma[1,1], length(df$group)+1)
      y2 <- rep(limits.2sigma[1,2], length(df$group)+1)
    } else
    {
      limits.2sigma <- limits.2sigma[plot.index$row, , drop = FALSE]
      y1 <- c(limits.2sigma[,1], limits.2sigma[length(df$group),1])
      y2 <- c(limits.2sigma[,2], limits.2sigma[length(df$group),2])
    }
    xp1 <- rep(x1, each=2)[-1]
    xp2 <- rep(x2, each=2)[-1]
    yp1 <- rep(y1, each=2)[-2*length(y1)]
    yp2 <- rep(y2, each=2)[-2*length(y2)]
    if(fill)
    { # fill the in-control area
      plot <- plot + 
        geom_polygon(data = data.frame(x = c(xp1,rev(xp2)), 
                                       y = c(yp1,rev(yp2))),
                     aes(x = .data[["x"]], 
                         y = .data[["y"]]), 
                     fill = adjustcolor(qcc.options("zones")$fill, alpha.f=0.2),
                     col = NA)
    } else
    {
      plot <- plot + 
        geom_step(data = data.frame(x = x1, y = y1),
                  aes(x = .data[["x"]], 
                      y = .data[["y"]]), 
                  lty = qcc.options("zones")$lty[3],
                  col = qcc.options("zones")$col[3])
      plot <- plot + 
        geom_step(data = data.frame(x = x2, y = y2),
                  aes(x = .data[["x"]], 
                      y = .data[["y"]]), 
                  lty = qcc.options("zones")$lty[3],
                  col = qcc.options("zones")$col[3])
    }
  }
  
  # draw center line
  plot <- plot + if(length(center) == 1) 
  {
    geom_hline(yintercept = center, 
               col = qcc.options("zones")$col[1]) 
  } else
  {
    geom_step(data = df, aes(x = .data[["group"]], 
                             y = .data[["center"]]),
              col = qcc.options("zones")$col[1])
  }

  if(chart.all & (!is.null(newstats)))
  {
    len.obj.stats <- length(stats)
    len.new.stats <- length(newstats)
    plot <- plot +
      geom_vline(xintercept = mean(groups[len.obj.stats+c(0,1)]), lty = 3) +
      annotate("text", 
               x = mean(c(groups[1], mean(groups[len.obj.stats+c(0,1)]))),
               y = max(ylim),
               label = "Calibration data", 
               hjust = 0.5, vjust = -0.5, size = 10 * 5/14) +
      annotate("text",
               x = mean(c(mean(groups[len.obj.stats+c(0,1)]), 
                          groups[len.obj.stats+len.new.stats])),
               y = max(ylim),
               label = "New data", 
               hjust = 0.5, vjust = -0.5, size = 10 * 5/14)
  }
  
  if(add.stats) 
  { 
    # write info at bottom
    tab_base <- ggplot() + 
      ggplot2::xlim(0,1) + ggplot2::ylim(0,1) + 
      theme_void() +
      theme(plot.background = element_rect(fill = qcc.options("bg.margin"),
                                           color = qcc.options("bg.margin")),
            plot.margin = margin(0.5, 0, 0.5, 0, unit = "lines"))
    
    text1 <- paste(paste0("Number of groups = ", length(statistics)),
                   paste0("Center = ", if(length(center) == 1) 
                     signif(center[1], digits) else "variable"),
                   paste0("StdDev = ", if(length(std.dev) == 1) 
                     signif(std.dev[1], digits) else "variable"), 
                   sep = "\n")
    text2 <- paste("",
                   paste0("LCL = ", if(length(unique(lcl)) == 1) 
                     signif(lcl[1], digits) else "variable"),
                   paste0("UCL = " ,if(length(unique(ucl)) == 1) 
                     signif(ucl[1], digits) else "variable"), 
                   sep = "\n")
    text3 <- paste("",
                   paste0("No. beyond limits = ", sum(violations == 1, na.rm=TRUE)),
                   paste0("No. violating runs = ", sum(violations > 1, na.rm=TRUE)),
                   sep = "\n")
    tab1 <- tab_base + 
      geom_text(aes(x = -Inf, y = Inf), label = text1, 
                hjust = 0, vjust = 1, size = 10 * 5/14)
    tab2 <- tab_base + 
      geom_text(aes(x = -Inf, y = Inf), label = text2, 
                hjust = 0, vjust = 1, size = 10 * 5/14)
    tab3 <- tab_base + 
      geom_text(aes(x = -Inf, y = Inf), label = text3, 
                hjust = 0, vjust = 1, size = 10 * 5/14)

    plot <- patchwork::wrap_plots(plot, tab1, tab2, tab3, 
                                  design = c("AAA\nBCD"),
                                  heights = c(0.85, 0.15), 
                                  widths = c(0.4, 0.3, 0.3))
  }
  
  return(plot)
}
