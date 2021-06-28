#' List of vectors to UpSet data converter
#'
#' Same to \code{\link[UpSetR]{fromList}}, but add rownames.
#'
#' @param input List of sets.
#' @return A data.frame contains binary data.
#' @export
from_list_to_upset <- function(input) {
  data <- UpSetR::fromList(input)
  elements <- unique(unlist(input))
  rownames(data) <- elements
  data
}

#' Query element counts of intersections.
#'
#' @param bin Binary data produced by \code{from_list_to_upset}
#' @export
intersection_counts <- function(bin, ...) {
  sets <- as.character(list(...))
  sum(apply(bin, 1, function(x) {
    all(x[sets] == 1) && sum(x) == length(sets)
  }))
}

#' Convert UpSet binary data to tidy data.frame
#'
#' @param bin UpSet binary data produced from \code{export}
#' @return A data.frame
#' @export
upset_dataframe <- function(bin) {
  bin %>%
    tibble::as_tibble(rownames = "gene") %>%
    tidyr::gather(-gene, key = "cluster", value = "is_diff_expr") %>%
    dplyr::filter(is_diff_expr == 1) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarize(clusters = list(cluster))
}

#' Plot UpSet diagram from tidy UpSet data.frame
#'
#' @param df A data.frame produced from \code{upset_dataframe}
#' @param n_intersections How many intersections to show.
#' @param x_expansion Expansion for x-axis, see \code{\link[ggplot2]{expansion}}
#' @param y_expansion Expansion for y-axis, see \code{\link[ggplot2]{expansion}}
#' @param plot_margin margin around entire plot (unit with the sizes of the top,
#'  right, bottom, and left margins)
#' @param ... pass to \code{\link[ggupset]{scale_x_upset}}
#' @export
upset_plot <- function(
  df,
  n_intersections = 15,
  x_expansion = ggplot2::expansion(c(0, 0), c(0.8, 0.8)),
  y_expansion = ggplot2::expansion(c(0, 0.1), c(0, 0)),
  plot_margin = ggplot2::margin(0.5, 0.5, 0.5, 1.5, unit = "cm"),
  ...
) {
  main_plot <- df %>%
    ggplot2::ggplot(ggplot2::aes(x = clusters)) +
    ggplot2::geom_bar() +
    ggplot2::geom_text(
      stat = "count", size = 2, vjust = -1,
      mapping = ggplot2::aes(label = ggplot2::after_stat(count))) +
    ggupset::scale_x_upset(
      expand = x_expansion,
      n_intersections = n_intersections,
      ...) +
    ggplot2::scale_y_continuous(expand = y_expansion) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.margin = plot_margin)

  main_plot
}

#' Plot up/down-regulated genes in one diagram
#'
#' @param up_list A list of up-regulated gene sets
#' @param down_list A list of down-regulated gene sets
#' @param ... pass to \code{\link[ggupset]{scale_x_upset}}
#' @inheritParams upset_plot
#' @export
upset_updown_regulated <- function(
  up_list, down_list,
  n_intersections = 15,
  x_expansion = ggplot2::expansion(c(0, 0), c(0.8, 0.8)),
  y_expansion = ggplot2::expansion(c(0, 0.1), c(0, 0)),
  plot_margin = ggplot2::margin(0.5, 0.5, 0.5, 1.5, unit = "cm"),
  ...
) {

  up_data <- from_list_to_upset(up_list) %>% upset_dataframe()
  down_data <- from_list_to_upset(down_list) %>% upset_dataframe()

  up_data$regulated <- factor("up", c("up", "down"))
  down_data$regulated <- factor("down", c("up", "down"))

  df_to_plot <- dplyr::bind_rows(up_data, down_data)

  main_plot <- df_to_plot %>%
    ggplot2::ggplot(ggplot2::aes(x = clusters, fill = regulated)) +
    ggplot2::geom_bar(
      width = 0.9,
      position = ggplot2::position_dodge(width = 0.9)) +
    ggplot2::geom_text(
      stat = "count", size = 2, vjust = -1,
      position = ggplot2::position_dodge(width = 0.9),
      mapping = ggplot2::aes(label = ggplot2::after_stat(count))) +
    ggupset::scale_x_upset(
      expand = x_expansion,
      n_intersections = n_intersections,
      ...) +
    ggplot2::scale_y_continuous(expand =  y_expansion) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.margin = plot_margin)

  main_plot
}