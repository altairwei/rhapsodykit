#ggplot2::theme_set(cowplot::theme_cowplot())

gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#' Plot heatmap across samples.
#'
#' @param x Data object to plot.
#' @param genes Genes to be queried.
#' @param samples Samples to be queried.
#' @param clusters Cell subpopulation to be queried.
#' @param groups Group infromation.
#' @inheritParams calculate_pseudo_bulk
#' @export
heatmap_cross_sample <- function(
  x, genes, samples = NULL, clusters = NULL, groups = NULL,
  type = c("counts", "logcounts", "cpm", "vstresiduals"),
  ...) {
  UseMethod("heatmap_cross_sample")
}

#' @describeIn heatmap_cross_sample Plot heatmap for Seurat object.
#' @method heatmap_cross_sample Seurat
#' @export
heatmap_cross_sample.Seurat <- function(
  x, genes, samples = NULL, clusters = NULL, groups = NULL,
  type = c("counts", "logcounts", "cpm", "vstresiduals"),
  ...) {
  type <- match.arg(type)
  sce <- prepare_muscat_sce(x)
  heatmap_cross_sample(
    sce,
    genes = genes,
    samples = samples,
    clusters = clusters,
    type = type,
    ...
  )
}

#' @describeIn heatmap_cross_sample Plot heatmap for SingleCellExperiment object.
#' @method heatmap_cross_sample SingleCellExperiment
#' @export
heatmap_cross_sample.SingleCellExperiment <- function(
  x, genes, samples = NULL, clusters = NULL, groups = NULL,
  type = c("counts", "logcounts", "cpm", "vstresiduals"),
  ...) {
  type <- match.arg(type)
  sce <- calculate_pseudo_bulk(x, type)
  pb_list <- abind::abind(as.list(sce@assays@data), along = 3)
  heatmap_cross_sample(
    pb_list,
    genes = genes,
    samples = samples,
    clusters = clusters,
    type = type,
    ...
  )
}

#' @describeIn heatmap_cross_sample Plot heatmap for a list of matrices.
#' @method heatmap_cross_sample list
#' @export
heatmap_cross_sample.list <- function(
  x, genes, samples = NULL, clusters = NULL, groups = NULL,
  type = c("counts", "logcounts", "cpm", "vstresiduals"),
  ...) {
  type <- match.arg(type)
  stopifnot(
    !is.null(names(x)),
    all(sapply(x, function(m) is.matrix(m))),
    !is.null(genes)
  )

  # Validate subset keys

  if (!is.null(clusters)) {
    clusters_valid <- clusters %in% names(x)
    if (!all(clusters_valid)) {
      warning("Not all queried cluster exist.")
      clusters <- clusters[clusters_valid]
    }
  } else {
    clusters <- names(x)
  }

  gene_names_list <- lapply(x, function(m) rownames(m))
  genes_all_identical <- all(sapply(
    gene_names_list, FUN = identical, gene_names_list[[1]]))
  if (!genes_all_identical)
    stop("Matrix must have same rownames.")

  sample_names_list <- lapply(x, function(m) colnames(m))
  sample_all_identical <- all(sapply(
    sample_names_list, FUN = identical, sample_names_list[[1]]))
  if (!sample_all_identical)
    stop("Matrix must have same colnames.")

  if (!is.null(samples)) {
    samples_valid <- samples %in% colnames(x[[1]])
    if (!all(samples_valid)) {
      warning("Not all queried samples exist.")
      samples <- samples[samples_valid]
      if (length(samples) == 0)
        stop("No samples available.")
    }
  } else {
    samples <- colnames(x[[1]])
  }

  genes_valid <- genes %in% gene_names_list[[1]]
  if (!all(genes_valid)) {
    warning("Not all queried genes exist.")
    genes <- genes[genes_valid]
    if (length(genes) == 0)
      stop("No genes available.")
  }

  # Subset raw data
  x <- x[clusters]
  x <- lapply(x, function(m) m[genes, samples, drop = FALSE])

  clu_fct <- factor(rep(names(x), each = length(genes)))

  mat_combined <- do.call(rbind, x)

  ComplexHeatmap::Heatmap(
    mat_combined,
    name = type,
    row_split = clu_fct,
    ...
  )
}

#' @describeIn heatmap_cross_sample Plot heatmap for array.
#' @method heatmap_cross_sample array
#' @export
heatmap_cross_sample.array <- function(
  x, genes, samples = NULL, clusters = NULL, groups = NULL,
  type = c("counts", "logcounts", "cpm", "vstresiduals"),
  ...) {
  type <- match.arg(type)
  stopifnot(
    length(dim(x)) == 3,
    length(dimnames(x)) == 3
  )

  if (any(sapply(dimnames(x), is.null)))
    stop("All dims must have dimnames.")

  if (!is.null(clusters)) {
    clusters_valid <- clusters %in% dimnames(x)[[3L]]
    if (!all(clusters_valid)) {
      warning("Not all queried cluster exist.")
      clusters <- clusters[clusters_valid]
    }
  } else {
    clusters <- dimnames(x)[[3L]]
  }

  if (!is.null(samples)) {
    samples_valid <- samples %in% colnames(x)
    if (!all(samples_valid)) {
      warning("Not all queried samples exist.")
      samples <- samples[samples_valid]
      if (length(samples) == 0)
        stop("No samples available.")
    }
  } else {
    samples <- colnames(x)
  }

  genes_valid <- genes %in% rownames(x)
  if (!all(genes_valid)) {
    warning("Not all queried genes exist.")
    genes <- genes[genes_valid]
    if (length(genes) == 0)
      stop("No genes available.")
  }

  # Subset raw data
  x <- x[genes, samples, clusters, drop = FALSE]
  clu_fct <- factor(rep(dimnames(x)[[3L]], each = length(genes)))

  x_list <- purrr::array_branch(x, margin = 3)
  mat_combined <- do.call(rbind, x_list)

  ComplexHeatmap::Heatmap(
    mat_combined,
    name = type,
    row_split = clu_fct,
    ...
  )
}

#' Volcano plot for differential state analysis.
#'
#' @param x A data.frame contains following column:
#' \describe{
#'   \item{\code{gene}}{gene names}
#'   \item{\code{logFC}}{log2 fold change}
#'   \item{\code{p_adj.loc}}{adjusted p-value}
#'   \item{\code{cluster_id}}{data source}
#' }
#' @param logfc_cut Where to add cutoff-line for logFC.
#' @param pval_cut Where to add cutoff-line for P value.
#' @param y_var Which column should be used as \code{y}.
#' @param title Plot title.
#' @export
volcano_diff_state <- function(
  x,
  logfc_cut = 1,
  pval_cut = 0.05,
  y_var = "p_adj.loc",
  title = NULL
) {
  # setting for color
  x$color_transparent <- ifelse(
    (x[[y_var]] < pval_cut & x$logFC > logfc_cut), "red",
    ifelse((x[[y_var]] < pval_cut & x$logFC < -logfc_cut), "blue", "grey")
  )
  # setting for size
  size <- ifelse((x[[y_var]] < pval_cut & abs(x$logFC) > logfc_cut), 1, 0.5)

  x$y <- -log10(x[[y_var]])

  # Construct the plot object
  p1 <- ggplot2::ggplot(x, ggplot2::aes(logFC, y)) +
    ggplot2::geom_point(
      size = size, colour = x$color_transparent) +
    ggplot2::labs(
      x = bquote(~Log[2]~"(fold change)"),
      y = bquote(~-Log[10]~italic("P-adjusted")), title = "") +
    ggplot2::geom_vline(xintercept = c(-logfc_cut, logfc_cut), color = "grey40",
              linetype = "longdash", lwd = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(pval_cut), color = "grey40",
              linetype = "longdash", lwd = 0.5) +

    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::facet_wrap(~cluster_id)

  if (!is.null(title)) {
    p1 <- p1 + ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }

  p1
}

#' Plot number of genes per subpopulation per contrast.
#'
#' @param results Results from \code{\link{pseudobulk_diff_state}}
#' @param simplify If \code{TRUE}, subpopulation without DS genes
#'  won't be shown.
#' @export
genecount_diff_state <- function(
  results, simplify = FALSE) {
  p <- results %>%
    diff_state_format() %>%
    ggplot2::ggplot(ggplot2::aes(cluster_id, fill = contrast)) +
    ggplot2::geom_bar(
      position = ggplot2::position_dodge2(preserve = "single", padding = 0)) +
    ggplot2::ylab("Number of DS genes") +
    ggplot2::xlab("Cell Subpopulation") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!simplify)
    p <- p + ggplot2::scale_x_discrete(limits = names(results$data))

  p
}

#' Plot Abundance Comparison Across Sample Clusters
#'
#' @export
barplot_cluster_abundance <- function(x, ...) {
  UseMethod("barplot_cluster_abundance")
}

#' Calculate the percentage of cell types in samples
#'
#' @param clusters A vector or factor which indicates cell types.
#' @param samples A vector which indicates sample id.
#' @param groups A vector which indicates group id.
#' @return A data frame
#' @export
calculateClusterProportion <- function(clusters, samples, groups) {
  cluster_prop <- prop.table(table(clusters, samples), 2)

  df <- as.data.frame(cluster_prop, stringsAsFactors = FALSE)
  names(df) <- c("cluster_id", "sample_id", "frequency")

  if (is.factor(clusters)) {
    df$cluster_id <- factor(df$cluster_id, levels = levels(clusters))
  }

  df$group_id <- groups[match(df$sample_id, samples)]

  df <- dplyr::select(df, cluster_id, sample_id, group_id, frequency)

  df
}

#' @describeIn calculateClusterProportion Calculate the size of cell
#' types in samples
#'
#' @export
calculateClusterSize <- function(clusters, samples, groups) {
  cluster_size <- table(clusters, samples)

  df <- as.data.frame(cluster_size, stringsAsFactors = FALSE)
  names(df) <- c("cluster_id", "sample_id", "size")

  if (is.factor(clusters)) {
    df$cluster_id <- factor(df$cluster_id, levels = levels(clusters))
  }

  df$group_id <- groups[match(df$sample_id, samples)]

  df <- dplyr::select(df, cluster_id, sample_id, group_id, size)

  df
}

#' @describeIn barplot_cluster_abundance Plot abundance comparison across
#' sample clusters from a formatted data.frame
#' @method barplot_cluster_abundance data.frame
#' @param df Data frame to plot, must contains columns:
#' \describe{
#'   \item{\code{frequency}}{the proportion of cells in a cluster in a sample}
#'   \item{\code{sample_id}}{IDs of samples}
#'   \item{\code{cluster_id}}{IDs of clusters}
#'   \item{\code{group_id}}{IDs of sample groups}
#' }
#' @param position Specify the way to display the bar chart
#' @return A ggplot2 object
#' @export
barplot_cluster_abundance.data.frame <- function(
  df, position = c("stack", "dodge")
) {
  position <- match.arg(position)

  if (position == "stack") {
    p <- df %>%
      ggplot2::ggplot(ggplot2::aes(
        x = sample_id, y = frequency, fill = cluster_id)) +
      ggplot2::facet_wrap(~group_id, ncol = 1, scales = "free_y") +
      ggplot2::geom_bar(
        stat = "identity", col = "white",
        width = 1, size = 0.2, position = position) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::theme(
        aspect.ratio = NULL,
        panel.grid = ggplot2::element_blank(),
        panel.spacing = grid::unit(1, "mm"),
        strip.text = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank()
      )
  } else if (position == "dodge") {
    p <- df %>%
      ggplot2::ggplot(ggplot2::aes(
        x = sample_id, y = frequency, fill = group_id)) +
      ggplot2::facet_wrap(~cluster_id) +
      ggplot2::geom_bar(
        stat = "identity", col = "white",
        width = 1, size = 0.2,
        position = ggplot2::position_dodge()) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::theme(
        aspect.ratio = NULL,
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank()
      )
  }

  p
}

#' @describeIn barplot_cluster_abundance Calculate the proportion of cells from
#' object SingleCellExperiment
#' @inheritParams calculate_pseudo_bulk
#' @method barplot_cluster_abundance SingleCellExperiment
#' @export
barplot_cluster_abundance.SingleCellExperiment <- function(sce, ...) {
  df <- calculateClusterProportion(
    sce$cluster_id, sce$sample_id, sce$group_id)

  barplot_cluster_abundance(df, ...)
}

#' @describeIn barplot_cluster_abundance Calculate the proportion of cells from
#' Seurat object.
#' @method barplot_cluster_abundance Seurat
#'
#' @param srt Seurat object. The \code{sample} and \code{group} infromation
#' must exist in \code{seurat_object@meta.data} slot.
#'
#' @export
barplot_cluster_abundance.Seurat <- function(srt, ...) {
  stopifnot(
    inherits(srt, "Seurat"),
    !is.null(srt@meta.data$sample),
    !is.null(srt@meta.data$group)
  )

  srt_df <- SeuratObject::FetchData(srt, c("ident", "sample", "group"))
  df <- calculateClusterProportion(
    srt_df$ident, srt_df$sample, srt_df$group)

  barplot_cluster_abundance(df, ...)
}