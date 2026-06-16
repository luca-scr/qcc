theme_qcc <- function(...) {
  bg_margin <- qcc.options("bg.margin")
  bg_figure <- qcc.options("bg.figure")

  theme_light() +
    theme(
      plot.background = element_rect(fill = bg_margin, color = bg_margin),
      panel.background = element_rect(fill = bg_figure),
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "none",
      plot.margin = margin(5, 30, 5, 5),
      axis.text.y = element_text(
        angle = 90,
        margin = margin(l = 5, r = 5),
        hjust = 0.5,
        vjust = 0.5
      )
    ) +
    theme(...)
}
