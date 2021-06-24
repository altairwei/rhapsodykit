#' Perform enrichment analysis by using clusterProfiler.
#'
#' @param genes A vector of gene id.
#' @param go_data A data.frame contains four columns:
#' \enumerate{
#'   \item First column must be named as \code{Gene_ID}
#'   \item Second column must be named as \code{GO_ID}
#'   \item Third column must be named as \code{GO_Name}
#'   \item Fourth column must be named as \code{GO_Level}
#' }
#' @export
enrichment_analysis <- function(genes, go_data) {
  # Discard genes without annotation
  go_data <- go_data %>%
    dplyr::filter(GO_ID != "")

  # Construct TERM2GENE and TERM2NAME
  # TERM2GENE is a data.frame with first column of term ID and second column of
  # corresponding mapped gene and TERM2NAME is a data.frame with first column of
  # term ID and second column of corresponding term name
  go2gene <- go_data %>%
    dplyr::select(GO_ID, Gene_ID)
  go2name <- go_data %>%
    dplyr::select(GO_ID, GO_Name)

  # Perform GO enrichment analysis
  enrich_results <- clusterProfiler::enricher(
    genes, TERM2GENE = go2gene, TERM2NAME = go2name)

  # Add GO_Level to results
  go_level_table <- go_data %>%
    dplyr::group_by(GO_ID) %>%
    dplyr::summarize(Level = unique(GO_Level))
  enrich_results_with_level <- enrich_results %>%
    tibble::as_tibble() %>%
    dplyr::left_join(go_level_table, by = c("ID" = "GO_ID")) %>%
    dplyr::select(ID, Description, Level, dplyr::everything())
  enrich_results_with_level <- enrich_results_with_level %>%
    # RichFactor means DEGs of a GO term divided by All Genes of this GO term.
    dplyr::mutate(
      RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) %>%
    dplyr::select(ID:GeneRatio, BgRatio, RichFactor, dplyr::everything())

  enrich_results_with_level
}

#' Make barplot of GO enrichment analysis.
#'
#' @param resdf Results returned by \code{enrichment_analysis}
#' @export
enrich_barplot <- function(
  resdf, x = "GeneRatio", show_category = 10, title = "") {
  # Filter empty terms
  resdf <- resdf[!is.na(resdf$Description), ]
  resdf <- resdf[resdf$Count != 0, ]

  # Calculate ratios
  resdf$GeneRatio <- parse_ratio(resdf$GeneRatio)
  resdf$BgRatio <- parse_ratio(resdf$BgRatio)

  # Get top N GO terms
  if (show_category <= nrow(resdf)) {
      resdf <- resdf[1:show_category, ]
  }

  # Make top Go terms get higher y values, because of reverse levels
  resdf$Description <- factor(resdf$Description,
                          levels = rev(unique(resdf$Description)))

  # Set plot arguments
  colorBy <- "p.adjust"
  size <- "Count"

  # Set x to user specified col and then
  resdf <- dplyr::mutate(resdf, x = eval(parse(text = x)))

  # Re-order Go term description by x values,
  idx <- order(resdf[["x"]], decreasing = TRUE)
  resdf$Description <- factor(
    resdf$Description, levels = rev(unique(resdf$Description[idx])))

  # Make the plot
  offset <- max(resdf[[x]]) * 0.01
  resdf[["p.adjust.log"]] <- -log10(resdf[["p.adjust"]])

  ggplot2::ggplot(resdf,
      ggplot2::aes_string(x, "Description")) +
      ggplot2::geom_col(
        ggplot2::aes_string(alpha = "p.adjust.log"), fill = "deepskyblue") +
      ggplot2::scale_color_viridis_c() +
      ggplot2::scale_size_continuous(range = c(2, 10)) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_alpha_continuous(name = "-Log10(p.adjust)") +
      ggplot2::ylab(NULL) +
      ggplot2::geom_text(
        ggplot2::aes_string(label = "Description"),
        x = offset, hjust = 0) +
      ggplot2::ggtitle(title) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank())
}

#' Plot barplot of each GO level.
#'
#' @inheritParams enrich_barplot
#' @export
enrich_barplot_by_level <- function(
  resdf, x = "GeneRatio", show_category = 20) {
  p_list <- list()
  mf <- resdf %>%
      dplyr::filter(Level == "molecular_function")
  if (nrow(mf) > 0) {
    p_list[[length(p_list) + 1]] <- mf %>% go_barplot(x = x,
      show_category = show_category, title = "Molecular Function")
  }

  bp <- resdf %>%
    dplyr::filter(Level == "biological_process")
  if (nrow(bp) > 0) {
    p_list[[length(p_list) + 1]] <- bp %>% go_barplot(x = x,
      show_category = show_category, title = "Biological Process")
  }

  cc <- resdf %>%
      dplyr::filter(Level == "cellular_component")
  if (nrow(cc) > 0) {
    p_list[[length(p_list) + 1]] <- cc %>% go_barplot(x = x,
      show_category = show_category, title = "Cellular Component")
  }

  p <- patchwork::wrap_plots(p_list)

  p
}

parse_ratio <- function(ratio) {
  ratio <- sub("^\\s*", "", as.character(ratio))
  ratio <- sub("\\s*$", "", ratio)
  numerator <- as.numeric(sub("/\\d+$", "", ratio))
  denominator <- as.numeric(sub("^\\d+/", "", ratio))
  return(numerator / denominator)
}