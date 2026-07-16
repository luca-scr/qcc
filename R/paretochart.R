#' Pareto chart
#'
#' Computes a table of statistics and plot a Pareto chart.
#'
#' A Pareto chart is a barplot where the categories are ordered in non
#' increasing order, and a line is also added to show the cumulative sum.
#'
#' @aliases paretoChart print.paretoChart plot.paretoChart blues.colors
#' @export paretoChart
#' @param data a vector of values. `names(data)` are used for labelling
#' the bars.
#' @param plot a logical specifying if the chart should be provided
#' (`TRUE`, default).
#' @param x an object of class `'paretoChart'` returned by a call to
#' `paretoChart()` function.
#' @param title a character string specifying the main title. Set `title =
#' NULL` to remove the title.
#' @param xlab a string specifying the label for the x-axis.
#' @param ylab a string specifying the label for the y-axis.
#' @param ylab2 a string specifying the label for the second y-axis on the
#' right side.
#' @param ylim a numeric vector specifying the limits for the y-axis.
#' @param col a value for the color, a vector of colors, or a palette for the
#' bars. See the help for [colors()] and [palette()].
#' @param ... catch other optional arguments.
#' @return Returns an object of class `'paretoChart'` containing the
#' descriptive statistics used to draw the Pareto chart. This object has
#' associated a `print` and `plot` method.
#' @author Luca Scrucca
#' @seealso [barplot()]
#' @references Mason, R.L. and Young, J.C. (2002) *Multivariate
#' Statistical Process Control with Industrial Applications*, SIAM.
#'
#' Montgomery, D.C. (2013) *Introduction to Statistical Quality Control*,
#' 7th ed. New York: John Wiley & Sons.
#'
#' Ryan, T. P. (2011), *Statistical Methods for Quality Improvement*, 3rd
#' ed. New York: John Wiley & Sons, Inc.
#'
#' Scrucca, L. (2004). qcc: an R package for quality control charting and
#' statistical process control. *R News* 4/1, 11-17.
#'
#' Wetherill, G.B. and Brown, D.W. (1991) *Statistical Process Control*.
#' New York: Chapman & Hall.
#' @keywords hplot
#' @examples
#'
#' defect  = c(80, 27, 66, 94, 33)
#' names(defect)  = c("price code", "schedule date", "supplier code", "contact num.", "part num.")
#' pc = paretoChart(defect, ylab = "Error frequency")
#' pc
#' plot(pc)
#'
#' plot(paretoChart(defect, ylab = "Error frequency"), col=rainbow(length(defect)))
#'
paretoChart <- function(data, ...)
{ 
  call <- match.call(expand.dots = TRUE)
  data.name <- deparse(substitute(data))
  data <- as.table(data)
  if(length(dim(data))>1) 
    stop("only one-dimensional object (table, vector, etc.) may be provided")
  # 
  data <- sort(data, decreasing = TRUE, na.last = TRUE)
  csum <- cumsum(data)
  tab <- cbind(data, csum, 
               data/max(csum, na.rm = TRUE)*100, 
               csum/max(csum, na.rm = TRUE)*100) 
  colnames(tab) <- c("Frequency", "Cum.Freq.", 
                     "Percentage", "Cum.Percent.")
  names(dimnames(tab)) <- c(data.name, "")
  
  # create object of class 'paretoChart'
  object <- list(call = call, 
                 data.name = data.name, 
                 tab = tab)
  class(object) <- "paretoChart"
  
  return(object)
}

#' @rdname paretoChart
#' @export
#' @export print.paretoChart
print.paretoChart <- function(x, digits = getOption("digits") - 3, ...)
{
  object <- x   # Argh.  Really want to use 'object' anyway
  cat(cli::rule(left = crayon::bold("Pareto Chart"), 
                width = min(getOption("width"),50)), "\n")
  print(object$tab, digits = digits, ...)
}

#' @rdname paretoChart
#' @export
#' @export plot.paretoChart
plot.paretoChart <- function(x, 
                             title, xlab,
                             ylab = "Frequency", 
                             ylab2 = "Cumulative percentage", 
                             ylim, 
                             col = blues.colors(nlevels),
                             ...)
{
  
  if(missing(title)) title <- paste("Pareto Chart for", x$data.name)
  if(missing(xlab)) xlab <- ""
  if(missing(ylim)) ylim <- c(0, max(x$tab[,"Cum.Freq."], na.rm = TRUE))
  
  df <- data.frame(x = rownames(x$tab),
                   f = x$tab[,"Frequency"],
                   p = x$tab[,"Cum.Freq."])
  df$x <- factor(df$x, levels = unique(df$x)[order(df$f, decreasing = TRUE)])
  nlevels <- nlevels(df$x)
  df$idx <- seq(nlevels)
  
  plot <- ggplot(data = df) +
    geom_bar(aes(x = .data[["x"]], 
                 y = .data[["f"]]),
             stat = "identity", fill = col) +
    geom_point(aes(x = .data[["idx"]], 
                   y = .data[["p"]]), 
               size = 2) +
    geom_line(aes(x = .data[["idx"]], 
                  y = .data[["p"]])) +
    labs(title = title, y = ylab, x = xlab) +
    scale_y_continuous(limits = ylim,
                       sec.axis = sec_axis(~./(max(.)*.95),
                                           name = ylab2,
                                           labels = scales::label_percent())) +
    theme_light() + 
    theme(plot.background = element_rect(fill = qcc.options("bg.margin"),
                                         color = qcc.options("bg.margin")),
          panel.background = element_rect(fill = qcc.options("bg.figure")),
          plot.title = element_text(face = "bold", size = 11),
          plot.margin = margin(5, 5, 5, 5))
  
  return(plot)
}
