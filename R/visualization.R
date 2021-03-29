plot_subpopulation_heatmap <- function(
  aggregated, show_gene_name = TRUE) {
  # Convert data frame to matrix
  # aggregated <- aggregated %>%
  #   tidyr::pivot_wider(names_from = ident, values_from = value)
  # mat <- data.matrix(aggregated[, -1])
  # rownames(mat) <- aggregated[["gene"]]

  p <- ggplot2::ggplot(aggregated, ggplot2::aes(ident, gene, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_fill_gradient(
      low = "white", high = "red") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45, vjust = 0, size = 12, hjust = 0),
      axis.text.y = ggplot2::element_text(size = 12)) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "top")

  if (isFALSE(show_gene_name))
    p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank())

  p
}