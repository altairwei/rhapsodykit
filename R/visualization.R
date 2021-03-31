plot_subpopulation_heatmap <- function(
                                       aggregated, show_gene_name = TRUE) {
  stopifnot(is.matrix(aggregated))

  ComplexHeatmap::Heatmap(aggregated,
    row_names_side = "left",
    column_names_side = "top",
    show_row_names = show_gene_name
  )

  # p <- ggplot2::ggplot(aggregated, ggplot2::aes(ident, gene, fill = value)) +
  #   ggplot2::geom_tile() +
  #   ggplot2::scale_x_discrete(position = "top") +
  #   ggplot2::scale_fill_gradient(
  #     low = "white", high = "red") +
  #   ggplot2::theme_minimal() +
  #   ggplot2::theme(
  #     axis.text.x = ggplot2::element_text(
  #       angle = 45, vjust = 0, size = 12, hjust = 0),
  #     axis.text.y = ggplot2::element_text(size = 12)) +
  #   ggplot2::theme(
  #     axis.title.x = ggplot2::element_blank(),
  #     axis.title.y = ggplot2::element_blank(),
  #     panel.grid.major = ggplot2::element_blank(),
  #     panel.border = ggplot2::element_blank(),
  #     panel.background = ggplot2::element_blank(),
  #     axis.ticks = ggplot2::element_blank(),
  #     legend.position = "top")

  # if (isFALSE(show_gene_name))
  #   p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank())

  # p
}

# Copy from InteractiveComplexHeatmap
get_pos_from_brush <- function(brush, ratio = 1) {
  coords <- brush$coords_css
  if (is.null(coords)) {
    return(NULL)
  }
  height <- (brush$range$bottom - brush$range$top) / brush$img_css_ratio$y
  pos1 <- grid::unit(c(coords$xmin, height - coords$ymin), "bigpts")
  pos2 <- grid::unit(c(coords$xmax, height - coords$ymax), "bigpts")
  pos1 <- pos1 / ratio
  pos2 <- pos2 / ratio
  list(pos1, pos2)
}

# Copy from InteractiveComplexHeatmap
get_pos_from_click <- function(click, ratio = 1) {
  if (identical(c("x", "y"), names(click))) {
    pos1 <- grid::unit(c(click$x, click$y), "bigpts")
  } else {
    coords <- click$coords_css
    if (is.null(coords)) {
      return(NULL)
    }
    height <- (click$range$bottom - click$range$top) / click$img_css_ratio$y
    pos1 <- grid::unit(c(coords$x, height - coords$y), "bigpts")
  }
  pos1[1] <- pos1[1] / ratio
  pos1[2] <- pos1[2] / ratio
  pos1
}