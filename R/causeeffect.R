#-------------------------------------------------------------------#
#                                                                   #
#                  Cause-and-Effect Diagram                         #
#                                                                   #
#-------------------------------------------------------------------#



#' Cause and Effect Diagram
#'
#' Draw a basic cause and effect diagram.
#'
#'
#' @param cause a list of causes and branches providing descriptive labels (see
#' the example below).
#' @param effect a string label or the effect.
#' @param title a character string specifying the main title. Set `title =
#' NULL` to remove the title.
#' @param cex a vector of values for the graphical character expansion. The
#' values refer, in order, to branches, causes and effect.
#' @param font a vector of values for the font to use. The values refer, in
#' order, to branches, causes and effect.
#' @param ... catches further ignored arguments.
#' @author Luca Scrucca
#' @references Montgomery, D.C. (2013) *Introduction to Statistical
#' Quality Control*, 7th ed. New York: John Wiley & Sons.
#'
#' Wetherill, G.B. and Brown, D.W. (1991) *Statistical Process Control*.
#' New York: Chapman & Hall.
#' @keywords hplot
#' @export
#' @examples
#'
#' causeEffectDiagram(cause = list(Measurements = c("Micrometers", 
#'                                                  "Microscopes", 
#'                                                  "Inspectors"),
#'                                 Materials = c("Alloys", 
#'                                               "Lubricants", 
#'                                               "Suppliers"),
#'                                 Personnel = c("Shifts", 
#'                                               "Supervisors", 
#'                                               "Training", 
#'                                               "Operators"),
#'                                 Environment = c("Condensation", 
#'                                                 "Moisture"),
#'                                 Methods = c("Brake",
#'                                             "Engager", 
#'                                             "Angle"),
#'                                 Machines = c("Speed", 
#'                                              "Lathes", 
#'                                              "Bits", 
#'                                              "Sockets")),
#'                    effect = "Surface Flaws")
#'
causeEffectDiagram <- function(cause, effect, 
                               title = "Cause-and-Effect diagram",
                               cex = c(0.9,1,1.2), 
                               font = c(3,1,2),
                               ...)
{
  
  # running mean of successive pairs of obs
  mean2 <- function(x)
  { 
    m <- rep(NA, length(x)-1)
    for (i in 1:(length(x)-1))
      m[i] <- mean(x[c(i,i+1)])
    return(m)
  }
  
  branch_count <- length(cause)
  upper_branch_count <- ceiling(branch_count/2)
  lower_branch_count <- branch_count - upper_branch_count
  max_subcauses_per_branch <- max(sapply(cause, length))

  plot <- ggplot() +
    xlim(0, 100) + ylim(0, 100) +
    labs(title = title) +
    theme_void() +
    theme(plot.background = element_rect(fill = qcc.options("bg.margin"),
                                         color = qcc.options("bg.margin")),
          panel.background = element_rect(fill = qcc.options("bg.figure")),
          plot.title = element_text(face = "bold", margin = margin(b = 10)),
          plot.margin = margin(10, 10, 10, 10))
  
  size <- cex*12/ggplot2::.pt
  inches_per_unit <- 6 / 100  # ~ 6 inches across 100 units
  effect_label_width <- strwidth(effect, units = "inches") /inches_per_unit * 1.1
  max_subcause_label_width <- max(unlist(sapply(cause, strwidth, units = "inches"))) /inches_per_unit # FIX: unused variable
  max_label_height <- max(strheight(effect, units="inches"),
            unlist(sapply(cause, strheight, units="inches"))) /inches_per_unit # FIX: unused variable

  plot <- plot +  
    # add main spine
    geom_segment(aes(x = 0, y = 50, xend = 100-effect_label_width-2, yend = 50),
                 arrow = arrow(length = unit(0.02, "npc")), 
                 linewidth = 1) +
    # add effect label
    annotate("text", x = 100-effect_label_width, y = 50,  
             label = effect, hjust = 0, 
             fontface = font[3], size = size[3])

  # draw branches and cause labels
  branch_slots <- max(upper_branch_count, lower_branch_count)
  branch_spacing <- (100-effect_label_width)/(branch_slots+1)
  ac <- branch_spacing*(0:branch_slots)
  branch_anchors <- mean2(ac)
  for (i in seq_len(upper_branch_count))
  {
    plot <- plot +
      annotate("segment", x = branch_anchors[i], y = 95, xend = ac[i+1], yend = 50) +
      annotate("text", x = branch_anchors[i], y = 97, label = names(cause)[i], 
               vjust = 0, fontface = font[2], size = size[2])
  }

  for (i in seq_len(lower_branch_count))
  {
    plot <- plot + 
      annotate("segment", x = branch_anchors[i], y = 5, xend = ac[i+1], yend = 50) +
      annotate("text", x = branch_anchors[i], y = 3, label = names(cause)[upper_branch_count+i], 
               vjust = 1, fontface = font[2], size = size[2])
  }
  # draw labels for upper branches
  for (j in 1:upper_branch_count)
  { 
    b <- (50-95)/(ac[j+1]-branch_anchors[j])
    intercept <- 95-b*branch_anchors[j]
    y <- rev(50+cumsum((95-50)/(max_subcauses_per_branch+1))*(1:(max_subcauses_per_branch)))
    x <- (y-intercept)/b
    for (i in 1:length(y))
    { 
      label <- cause[[j]][i]
      if (!is.na(label))
      {
        plot <- plot + 
          annotate("text", x = x[i], y = y[i], label = label, 
                   hjust = -0.1, fontface = font[1], size = size[1])
      }
    }
  }
  # draw labels for lower branches
  for (j in seq_len(lower_branch_count))
  { 
    b <- (50-5)/(ac[j+1]-branch_anchors[j])
    intercept <- 5-b*branch_anchors[j]
    y <- cumsum((95-50)/(max_subcauses_per_branch+1))*(1:(max_subcauses_per_branch))
    x <- (y-intercept)/b
    for (i in 1:length(y))
    {
      label <- cause[[upper_branch_count+j]][i]
      if (!is.na(label))
      { 
        plot <- plot + 
          annotate("text", x = x[i], y = y[i], label = label, 
                   hjust = -0.1, fontface = font[1], size = size[1])
      }
    }
  }
  
  return(plot)
}
