plot_subpopulation_heatmap <- function(
  aggregated, show_gene_name = TRUE) {
  stopifnot(is.matrix(aggregated))

  ComplexHeatmap::Heatmap(aggregated,
    row_names_side = "left",
    column_names_side = "top",
    show_row_names = show_gene_name
  )
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

plot_placeholder <- function(text) {
  grid::grid.text(
    text, 0.5, 0.5, gp = grid::gpar(fontsize = 20))
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
#' @param results results from \code{\link{pseudobulk_diff_state}}
#' @inheritParams diff_state_filter
#' @export
genecount_diff_state <- function(results, fdr_limit = 0.05, logfc_limit = 1) {
  results %>%
    lapply(function(x) {
      x %>%
        diff_state_filter(fdr_limit, logfc_limit) %>%
        dplyr::bind_rows()
    }) %>%
    dplyr::bind_rows() %>%
    ggplot2::ggplot(ggplot2::aes(cluster_id, fill = contrast)) +
      ggplot2::geom_bar(
        position = ggplot2::position_dodge2(preserve = "single", padding = 0)) +
      ggplot2::ylab("Number of DS genes") +
      ggplot2::xlab("Cell Subpopulation") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
