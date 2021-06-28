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
#' @param ... Arguments passed to \code{\link[clusterProfiler]{enricher}}
#' @return A \code{enrichResult} instance
#' @export
enrichment_analysis <- function(genes, go_data, ...) {
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
    genes, TERM2GENE = go2gene, TERM2NAME = go2name, ...)

  enrich_results
}

#' Add GO level data to result of enrichment analysis.
#'
#' @param enr A \code{enrichResult} instance
#' @inheritParams enrichment_analysis
#' @return A data.frame
#' @export
enrich_add_go_level <- function(enr, go_data) {
  # Discard genes without annotation
  go_data <- go_data %>%
    dplyr::filter(GO_ID != "")

  # Add GO_Level to results
  go_level_table <- go_data %>%
    dplyr::group_by(GO_ID) %>%
    dplyr::summarize(Level = unique(GO_Level))

  enrich_results_with_level <- enr %>%
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

#' Make barplot of enrichment analysis.
#'
#' @param obj Results returned by \code{enrichment_analysis}
#' @param x Which column used as x-axis
#' @param show_category How many annotation terms to show.
#' @param title Title of plot.
#' @export
enrich_barplot <- function(
  obj, x = "GeneRatio", show_category = 10, title = "") {
  UseMethod("enrich_barplot")
}

#' @describeIn enrich_barplot Barplot for \code{enrichResult} instance
#' @method enrich_barplot enrichResult
#' @export
enrich_barplot.enrichResult <- function(
  enr, x = "GeneRatio", show_category = 10, title = "") {
  enrich_barplot(
    enr@result, x = x, show_category = show_category, title = title)
}

#' @describeIn enrich_barplot Barplot for \code{data.frame}
#' @method enrich_barplot data.frame
#' @export
enrich_barplot.data.frame <- function(
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
      cowplot::theme_cowplot() +
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

#' Make Gene-Concept Network plot for enrichment analysis
#'
#' @param enr A \code{enrichResult} instance.
#' @param fold_change A named numeric vector of fold Change.
#'  The names of \code{fold_change} should be consistent with
#'  ID of \code{enr}
#' @param show_category How many annotation terms to show.
#' @param circular Whether using circular layout.
#' @export
enrich_cnetplot <- function(
  enr,
  fold_change = NULL,
  show_category = 10,
  circular = FALSE
) {
  p <- enr %>%
    enrichplot::cnetplot(
      foldChange = fold_change,
      showCategory = show_category,
      circular = circular,
      node_label = "none",
      colorEdge = circular
    )

  p <- p + ggraph::geom_node_text(
    ggplot2::aes_(label = ~name),
    data = p$data[1:show_category, ],
    repel = TRUE,
    bg.color = "white"
  )

  p
}

#' Prepare significant logFC vector named with gene names.
#'
#' @param ds Results returned by \code{\link[muscat]{pbDS}} or
#'  \code{\link{pseudobulk_diff_state}}
#' @inheritParams diff_state_pull
#' @param FDR Upper limit of FDR.
#' @param logFC Lower limit of absolute value of logFC.
#' @return A logFC vector named gene names
#' @export
enrich_significant_lfc <- function(ds, contrasts, clusters, FDR, logFC) {
  df <- ds %>%
      diff_state_significant(FDR, logFC) %>%
      diff_state_pull(contrasts, clusters, c("gene", "logFC"))

  fc_list <- df[, "logFC"]
  names(fc_list) <- as.character(df[, "gene"])

  fc_list
}

#' Prepare ranked logFC vector named with gene names.
#'
#' @inheritParams enrich_significant_lfc
#' @return A ranked logFC vector named gene names
#' @export
enrich_ranked_lfc <- function(ds, contrasts, clusters) {
  df <-  ds %>%
    # We don't need to filter genes with given cutoff
    diff_state_pull(contrasts, clusters, c("gene", "logFC"))

  fc_list <- df[, "logFC"]
  names(fc_list) <- as.character(df[, "gene"])

  ranked_list <- sort(fc_list, decreasing = TRUE)

  ranked_list
}

#' Perform over representation analysis
#'
#' @inheritParams enrich_significant_lfc
#' @inheritParams enrichment_analysis
#' @inheritDotParams enrichment_analysis
#' @return A \code{enrichResult} instance
#' @export
enrich_perform_ora <- function(
  ds, go_data, contrasts, clusters, FDR, logFC, ...
) {
  gene_lfc <- enrich_significant_lfc(ds, contrasts, clusters, FDR, logFC)

  go_data <- go_data %>%
    dplyr::filter(GO_ID != "")

  go2gene <- go_data %>%
    dplyr::select(GO_ID, Gene_ID)
  go2name <- go_data %>%
    dplyr::select(GO_ID, GO_Name)

  ora <- clusterProfiler::enricher(
    names(gene_lfc), TERM2GENE = go2gene, TERM2NAME = go2name, ...)

  ora
}

#' Perform gene set enrichment analysis
#'
#' @inheritParams enrich_ranked_lfc
#' @inheritParams enrichment_analysis
#' @param ... pass to \code{\link[clusterProfiler]{GSEA}}
#' @return A \code{gseaResult} instance
#' @export
enrich_perform_gsea <- function(ds, go_data, contrasts, clusters, ...) {
  ranked_list <- enrich_ranked_lfc(ds, contrasts, clusters)

  go_data <- go_data %>%
    dplyr::filter(GO_ID != "")

  go2gene <- go_data %>%
    dplyr::select(GO_ID, Gene_ID)
  go2name <- go_data %>%
    dplyr::select(GO_ID, GO_Name)

  gsea <- clusterProfiler::GSEA(
    ranked_list, TERM2GENE = go2gene, TERM2NAME = go2name, ...)

  gsea
}

#' Perform GSEA on gene clusters.
#'
#' @param ranked_clusters a list of ranked logFC vector.
#' @param ... pass to \code{\link[clusterProfiler]{GSEA}}
#' @return A un-official \code{compareClusterResult} object
#' @export
enrich_compare_gsea <- function(ranked_clusters, ...) {
  gsea_list <- lapply(ranked_clusters, function(i) {
    x <- suppressMessages(clusterProfiler::GSEA(i, ...))
    if (class(x) == "gseaResult") {
      as.data.frame(x)
    }
  })

  clusters_levels <- names(ranked_clusters)

  df <- dplyr::bind_rows(gsea_list, .id = "Cluster")
  df[["Cluster"]] <- factor(df$Cluster, levels = clusters_levels)

  new("compareClusterResult",
    compareClusterResult = df,
    geneClusters = ranked_clusters,
    fun = "GSEA",
    .call = match.call(expand.dots = TRUE)
  )
}

#' ep_str_wrap internal string wrapping function
#' @param string the string to be wrapped
#' @param width the maximum number of characters before wrapping to a new line
#' @author Guangchuang Yu
#' @noRd
ep_str_wrap <- function(string, width) {
  x <- gregexpr(" ", string)
  vapply(seq_along(x),
    FUN = function(i) {
      y <- x[[i]]
      n <- nchar(string[i])
      len <- (c(y, n) - c(0, y))
      idx <- len > width
      j <- which(!idx)
      if (length(j) && max(j) == length(len)) {
        j <- j[-length(j)]
      }
      if (length(j)) {
        idx[j] <- len[j] + len[j + 1] > width
      }
      idx <- idx[-length(idx)]
      start <- c(1, y[idx] + 1)
      end <- c(y[idx] - 1, n)
      words <- substring(string[i], start, end)
      paste0(words, collapse = "\n")
    },
    FUN.VALUE = character(1)
  )
}

#' Compare dotplot for GSEA
#'
#' @param object input object
#' @param show_category number of enriched terms to display
#' @param font_size font size
#' @param label_format a numeric value sets wrap length, alternatively
#'  a custom function to format axis labels. by default wraps names
#'  longer that 30 characters
#' @param title plot title
#' @export
dotplot_for_gsea <- function(
  object,
  show_category = 5,
  font_size = 12,
  label_format = 30,
  title = ""
) {
  df <- as.data.frame(object)

  if (is.null(show_category)) {
    result <- df
  } else if (is.numeric(show_category)) {
    top_n <- function(res, show_category) {
      res %>%
        dplyr::group_split(Cluster) %>%
        purrr::map_dfr(function(df, N) {
          if (length(df$setSize) > N) {
            idx <- order(df$pvalue, decreasing = FALSE)[1:N]
            return(df[idx, ])
          } else {
            return(df)
          }
        }, N = show_category)
    }

    result <- top_n(df, show_category)
  } else {
    result <- subset(df, Description %in% show_category)
  }

  ## remove zero count
  result$Description <- as.character(result$Description) ## un-factor
  GOlevel <- result[, c("ID", "Description")] ## GO ID and Term
  GOlevel <- unique(GOlevel)

  result <- result[result$setSize != 0, ]
  result$Description <- factor(
    result$Description, levels = rev(GOlevel$Description))

  core_count <- function(core_enrichment) {
    core_enrichment %>%
      stringr::str_split("/") %>%
      sapply(length)
  }

  result <- dplyr::mutate(result, Count = core_count(core_enrichment))

  label_func <- function(str) {
    ep_str_wrap(str, label_format)
  }

  if(is.function(label_format)) {
      label_func <- label_format
  }

  p <- ggplot2::ggplot(result, ggplot2::aes_string(
      x = "Cluster", y = "Description", size = "Count")) +
    ggplot2::geom_point(ggplot2::aes_string(color = "p.adjust")) +
    ggplot2::scale_color_continuous(
      low = "red", high = "blue",
      guide = ggplot2::guide_colorbar(reverse=TRUE)) +
    ggplot2::ylab(NULL) +
    ggplot2::ggtitle(title) +
    DOSE::theme_dose(font_size) +
    ggplot2::scale_size_continuous(range = c(3, 8)) +
    ggplot2::scale_y_discrete(labels = label_func)

  p
}