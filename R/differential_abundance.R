#' Run DAseq algorithm
#'
#' This function warp Steps 1~3 of DAseq algorithm from
#' package \code{\link[DAseq]{getDAcells}}
#'
#' @inheritParams DAseq::updateDAcells
#' @inheritParams DAseq::getDAcells
#' @inheritParams DAseq::getDAregion
#' @param obj Seurat object
#' @param ctrl Control condition
#' @param stim Stimulus condition
#' @param reduction Which dimensionality reduction to use
#' @param npc How many dimensions to use as input features
#'
#' @export
findDiffAbundantCells <- function(
  obj, ctrl, stim, reduction = "pca",
  npc = 20, pred.thres = NULL,
  k.vector = seq(50, 500, 50),
  resolution = 0.05
) {
  # DAseq only allow pairwise comparison, so we subset Seurat object.
  obj <- subset(obj, subset = group %in% c(ctrl, stim))
  info <- Seurat::FetchData(obj, vars = c("sample", "group", "ident"))

  # Extract sample label for each group.
  sample_ctrl <- unique(info[
    info$group %in% ctrl, "sample"])
  sample_stim <- unique(info[
    info$group %in% stim, "sample"])

  # Input data
  emb_pca <- Seurat::Embeddings(obj[[reduction]])[, 1:npc]

  # Sample label for every cell
  cell_labels <- obj$sample

  da_cells <- DAseq::getDAcells(
    X = emb_pca,
    cell.labels = cell_labels,
    labels.1 = sample_ctrl,
    labels.2 = sample_stim,
    k.vector = k.vector,
    do.plot = FALSE
  )

  if (!is.null(pred.thres))
    da_cells <- DAseq::updateDAcells(
      X = da_cells, pred.thres = pred.thres)

  da_regions <- DAseq::getDAregion(
    X = emb_pca,
    da.cells = da_cells,
    cell.labels = cell_labels,
    labels.1 = sample_ctrl,
    labels.2 = sample_stim,
    resolution = resolution,
    do.plot = FALSE
  )

  list(
    embeddings = list(
      pca = emb_pca,
      tsne = Seurat::Embeddings(obj[["tsne"]]),
      umap = Seurat::Embeddings(obj[["umap"]])
    ),
    cells = da_cells,
    regions = da_regions,
    info = info,
    args = list(
      cell.labels = cell_labels,
      labels.1 = sample_ctrl,
      condition.1 = ctrl,
      labels.2 = sample_stim,
      condition.2 = stim
    )
  )
}

#' Combine differential abundant subpopulation
#' TODO: 如何处理相同的细胞呢？去除重复，但用 list 记录该细胞来源条件。
#' @inherit Seurat::FindClusters
findDACombinedClusters <- function(obj_list, resolution) {

}

#' Find markers for differential abundant subpopulation
#'
#' @param obj Seurat object
#' @param da DAseq results
#' @param method Which method should be used to find markers.
#' Method \code{COSG} is based on \code{\link[COSG]{cosg}}, \code{Seruat} is
#' based on \code{\link[DAseq]{SeuratMarkerFinder}} and \code{STG} is baed
#' on \code{\link[DAseq]{STGmarkerFinder}} which requires python environment.
#' @param top_n Return the top N markers from results. The order of markers
#' is determined by scores for \code{COSG} and by average log2FC then P value
#' for \code{Seurat} or \code{STG}.
#' @param return_raw Return raw results produced by specified \code{method}.
#' @param GPU Which GPU to use (GPU IDs), default using CPU. Note: this value
#' will be used to set CUDA_VISIBLE_DEVICES environment.
#' @param ... Additional arguments passed to marker finder.
#' @export
findDiffAbundantMarkers <- function(
  obj, da,
  method = c("COSG", "Seurat", "STG"),
  top_n = NULL,
  return_raw = FALSE,
  GPU = "",
  ...
) {
  method <- match.arg(method)
  results <- switch(method,
    COSG = {
      obj <- DAseq::addDAslot(obj,
        da.regions = da$regions,
        da.slot = "da", set.ident = TRUE)
      n.da <- length(unique(obj$da)) - 1
      da.regions.to.run <- c(1:n.da)
      obj <- subset(x = obj, idents = as.character(da.regions.to.run))
      markers <- COSG::cosg(
        obj, groups = "all",
        assay = "RNA", slot = "data",
        n_genes_user = top_n,
        ...
      )

      if (return_raw)
        markers
      else
        as.list(markers$names)
    },
    Seurat = {
      obj <- DAseq::addDAslot(obj, da.regions = da$regions, da.slot = "da")
      markers <- DAseq::SeuratMarkerFinder(
        obj, da.slot = "da", assay = "RNA", ...
      )

      if (return_raw)
        markers
      else
        lapply(markers, function(df) {
          df <- dplyr::arrange(df, dplyr::desc(avg_log2FC), p_val)
          rownames(df)
        })
    },
    STG = {
      markers <- DAseq::STGmarkerFinder(
        X = SeuratObject::GetAssayData(
          obj, slot = "data", assay = "RNA"),
        da.regions = da$regions,
        return.model = TRUE,
        GPU = GPU,
        python.use = Sys.which("python"),
        ...
      )

      if (return_raw)
        markers
      else
        lapply(markers$da.markers, function(df) {
          dplyr::arrange(df, dplyr::desc(avg_logFC), p_value) %>%
          dplyr::pull("gene")
        })
    }
  )

  if (!return_raw && !is.null(top_n)) {
    lapply(results, function(x) x[seq_len(top_n)])
  } else {
    results
  }
}

#' Calculate DA Score
#'
#' @param cell.labels Sample label for all cells.
#' @param cell.idx Location of cells belong to a given identity.
#' @param labels.1 label name(s) that represent condition 1
#' @param labels.2 label name(s) that represent condition 2
#' @export
onlyDAscore <- function(cell.labels, cell.idx, labels.1, labels.2){
  # Remove invalid conditions
  labels.1 <- labels.1[labels.1 %in% cell.labels]
  labels.2 <- labels.2[labels.2 %in% cell.labels]

  # Get cell labels belong to one identity
  idx.label <- cell.labels[cell.idx]

  # This is same way to calculate DA.score as DA subpopulation
  ratio.1 <- sum(idx.label %in% labels.1) / sum(cell.labels %in% labels.1)
  ratio.2 <- sum(idx.label %in% labels.2) / sum(cell.labels %in% labels.2)
  ratio.diff <- (ratio.2 - ratio.1) / (ratio.2 + ratio.1)

  return(ratio.diff)
}

#' Plot DA Cell Scores
#'
#' The prediction values are overlayed on the 2D embedding
#'
#' @param obj DAseq results
#' @param show_clusters Display cell clusters
#' @param label Label cell custers
#' @param label_box Show label in box
#' @param label_cols Label colors
#' @param label_legend Show label legend
#' @return ggplot object
#' @export
plotDACellScore <- function(
  obj, reduction,
  show_clusters = NULL,
  label = FALSE,
  label_box = FALSE,
  label_cols = NULL,
  label_legend = TRUE,
  circle = FALSE
) {
  # Prepare data
  embedding <- obj$embeddings[[reduction]][obj$cells$cell.idx,]
  data_to_plot <- data.frame(
    Dim1 = embedding[, 1],
    Dim2 = embedding[, 2],
    clusters = obj$info$ident[obj$cells$cell.idx],
    score = obj$cells$da.pred
  )

  if (!is.null(show_clusters)) {
    data_to_plot <- dplyr::filter(data_to_plot, clusters %in% show_clusters)
    if (!is.null(label_cols)) {
      names(label_cols) <- levels(data_to_plot$clusters)
      label_cols <- label_cols[show_clusters]
    }
  }

  p <- data_to_plot %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(
      mapping = ggplot2::aes(x = Dim1, y = Dim2, col = score),
      shape = 16, size = 0.5
    ) +
    ggplot2::scale_color_gradientn(colours = c("blue", "white", "red")) +
    #ggthemes::scale_color_gradient2_tableau(trans = "reverse") +
    cowplot::theme_cowplot() +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  if (label) {
    label_positions <- data_to_plot %>%
      dplyr::group_by(clusters) %>%
      dplyr::summarise(x = median(Dim1), y = median(Dim2))

    p <- p +
      ggnewscale::new_scale_color() +
      ggrepel::geom_label_repel(
        data = label_positions,
        mapping = ggplot2::aes(
          x = x, y = y,
          label = clusters,
          color = clusters
        ),
        fill = ggplot2::alpha(c("white"), 0.8),
        show.legend = label_legend
      )

    if (!is.null(label_cols)) {
      p <- p + ggplot2::scale_color_manual(values = label_cols)
    }
  }

  if (circle) {
    p <- p + ggplot2::stat_ellipse(
      mapping = ggplot2::aes(x = Dim1, y = Dim2, group = clusters)
    )
  }

  p
}

#' Plot DA Sample Labels
#'
#' @inheritParams plotDACellScore
#' @export
plotDASampleLabel <- function(obj, reduction) {
  embedding <- obj$embeddings[[reduction]]
  data_to_plot <- data.frame(
    Dim1 = embedding[, 1],
    Dim2 = embedding[, 2],
    Group = obj$info$group
  )

  p <- data_to_plot %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(
      mapping = ggplot2::aes(x = Dim1, y = Dim2, col = Group),
      shape = 16, size = 0.5
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 3))) +
    cowplot::theme_cowplot() +
    ggplot2::theme(legend.title = ggplot2::element_blank())
  p
}

#' Selected DA cells are highlighted in the 2D embedding
#'
#' @inheritParams plotDACellScore
#' @export
plotDASite <- function(obj, reduction) {
  embedding <- obj$embeddings[[reduction]][obj$cells$cell.idx, ]
  DAseq::plotDAsite(embedding,
    site.list = list(obj$cells$da.down, obj$cells$da.up),
    cols = c("blue", "red")
  )
}

#' DA subpopulation clustering result is shown in the 2D embedding
#'
#' @inheritParams plotDACellScore
#' @export
plotDACellLabel <- function(obj, reduction) {
  embedding <- obj$embeddings[[reduction]][obj$cells$cell.idx, ]
  X.da.label <- obj$regions$da.region.label
  X.da.order <- order(X.da.label, decreasing = FALSE)
  X.n.da <- length(unique(X.da.label)) - 1
  DAseq::plotCellLabel(
    X = embedding[X.da.order,],
    label = as.factor(X.da.label[X.da.order])
  ) +
    ggplot2::scale_color_manual(
      values = c("gray", scales::hue_pal()(X.n.da))
    )
}

#' Plot the Overlap Between Cell Type and DA Subpopulations
#'
#' @inheritParams plotDACellScore
#' @export
plotDAOverlap <- function(obj) {
  contrib_df <- table(
      da.region.label = obj$regions$da.region.label,
      cell.idents = obj$info$ident) %>%
    apply(2, function(x) x / sum(x)) %>%
    as.data.frame.table()

  contrib_df$da.score <- apply(contrib_df, 1, function(x) {
    onlyDAscore(
      cell.labels = obj$args$cell.labels, 
      labels.1 = obj$args$labels.1, 
      labels.2 = obj$args$labels.2, 
      cell.idx = which(
        x[["cell.idents"]] == obj$info$ident &
        x[["da.region.label"]] == obj$regions$da.region.label
      )
    )
  })

  contrib_df %>%
    dplyr::filter(Freq != 0) %>%
    ggplot2::ggplot(ggplot2::aes(
      cell.idents, da.region.label,
      size = Freq, color = da.score)) +
    ggplot2::geom_point() +
    ggplot2::scale_y_discrete(
      labels = function(x) ifelse(x == 0, 0, paste0("DA", x))) +
    ggplot2::scale_color_gradientn(
      colours = c("blue", "white", "red"), limits = c(-1, 1)) +
    ggplot2::scale_size_continuous(labels = scales::percent) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    NULL
}

#' Plot the Random Permutation of DA Measure
#'
#' @inheritParams plotDACellScore
#' @param size Point size
#' @export
plotDARandomPermutation <- function(obj, size = 0.5) {
  n.cells <- length(obj$args$cell.labels)
  X.random.pred <- unlist(obj$cells$rand.pred)
  X.pred <- obj$cells$da.pred

  X.rand.plot <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = data.frame(
        order = seq(1, n.cells, length.out = length(X.random.pred)),
        random = sort(X.random.pred)
      ),
      mapping = ggplot2::aes(order, random),
      col = "gray", size = size, alpha = 0.5
    ) +
    ggplot2::geom_point(
      data = data.frame(
        order = c(1:n.cells),
        da = sort(X.pred)
      ),
      mapping = ggplot2::aes(order, da),
      col = "black", size = size, alpha = 0.75
    ) +
    ggplot2::geom_hline(yintercept = min(X.random.pred), size = size) +
    ggplot2::geom_hline(yintercept = max(X.random.pred), size = size) +
    ggplot2::scale_y_continuous(n.breaks = 11, limits = c(-1, 1)) +
    cowplot::theme_cowplot() +
    ggplot2::xlab(ggplot2::element_blank()) +
    ggplot2::ylab("DA measure")

  X.rand.plot
}
