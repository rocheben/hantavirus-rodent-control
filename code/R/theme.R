## Shared plotting theme and palettes ----------------------------------------

seov_theme <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92"),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      legend.title     = ggplot2::element_text(size = base_size - 1),
      legend.text      = ggplot2::element_text(size = base_size - 2),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size)
    )
}

## Colours of the three commensal:wild abundance ratios (Fig. 2): a single
## sequential ramp, dark where commensal species dominate.
seov_ratio_colours <- c("5" = "#7A0177", "1" = "#DD3497", "0.5" = "#FA9FB5")
