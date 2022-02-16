# Set theme
theme_rhr <-
  ggplot2::theme_bw(
    base_family = "Helvetica",
    base_size = 10
  ) +
  ggplot2::theme(
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(vjust = 0.6),
    panel.spacing = unit(0.1, "lines")
  )
