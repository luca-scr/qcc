#-------------------------------------------------------------------#
#                                                                   #
#                     PARETO CHART                                  #
#                                                                   #
#-------------------------------------------------------------------#

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

print.paretoChart <- function(x, digits = getOption("digits") - 3, ...)
{
  object <- x   # Argh.  Really want to use 'object' anyway
  cat(cli::rule(left = crayon::bold("Pareto Chart"), 
                width = min(getOption("width"),50)), "\n")
  print(object$tab, digits = digits, ...)
}

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
    geom_bar(aes_string(x = "x", y = "f"),
             stat = "identity", fill = col) +
    geom_point(aes_string(x = "idx", y = "p"), size = 2) +
    geom_line(aes_string(x = "idx", y = "p")) +
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
