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
  k.vector = NULL,
  resolution = 0.05
) {
  # DAseq only allow pairwise comparison, so we subset Seurat object.
  obj <- subset(obj, subset = group %in% c(ctrl, stim))
  info <- Seurat::FetchData(obj, vars = c("sample", "group", "ident"))

  # Extract sample label for each group. Aka. replicates
  sample_ctrl <- unique(info[
    info$group == ctrl, "sample"])
  sample_stim <- unique(info[
    info$group == stim, "sample"])

  # Input data
  emb_pca <- Seurat::Embeddings(obj[[reduction]])[, 1:npc]

  # Sample label for every cell
  cell_labels <- obj$sample

  # k.vector here means different nearest-neighborhood scales
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
#'
#' This function assumes that each DAseq result is produced
#' from a different subset of the same integrated cell embedding
#' and has a common reference condition (aka. control).
#'
#' @inheritParams DAseq::getDAregion
#' @param obj_list A list of DAseq results
#' @param ref_label The control condition common to all DAseq results
#' @export
findDACombinedClusters <- function(
  obj_list, ref_label,
  resolution = 0.05, prune.SNN = 1 / 15,
  group.singletons = FALSE,
  min.cell = NULL,
  ...) {

  stopifnot(length(obj_list) > 1)
  if (!all(sapply(obj_list, function(x) x$args$condition.1) == ref_label))
    stop("All DAseq results must have the same reference label: ", ref_label)

  # Combine cell embedings
  # The order of cell.labels is equal to the order of cell embeding
  embed_ref <- obj_list[[1]]$embeddings$pca[
    with(obj_list[[1]]$args, cell.labels %in% labels.1), ]
  embed_stim_list <- lapply(obj_list, function(obj) {
    obj$embeddings$pca[
      with(obj$args, cell.labels %in% labels.2), ]
  })

  embed <- rbind(embed_ref, do.call(rbind, embed_stim_list))
  if (nrow(embed) != length(unique(rownames(embed))))
    stop("Duplicated cell names.")

  # Generate meta data
  df_list <- lapply(obj_list, function(obj) {
    with(obj, {
      df <- data.frame(
        cell.names = rownames(embeddings$pca),
        sample = args$cell.labels,
        da.up = FALSE,
        da.down = FALSE
      )

      # The value of da.up or da.down is the index of cell.labels
      df$da.up[cells$da.up] <- TRUE
      df$da.down[cells$da.down] <- TRUE

      df
    })
  })

  meta.data <- dplyr::bind_rows(df_list) %>%
    dplyr::group_by(cell.names) %>%
    dplyr::summarise(
      cell.names = unique(cell.names),
      sample = unique(sample),
      da.up = any(da.up),
      da.down = any(da.down)
    )

  meta.data <- tibble::column_to_rownames(meta.data, var = "cell.names")
  meta.data <- meta.data[rownames(embed), ]

  n.cells <- nrow(embed)
  n.dims <- ncol(embed)
  X.S <- SeuratObject::CreateSeuratObject(
    counts = t(embed), meta.data = meta.data)
  # We assign it to 'pca' even though it may not be, it doesn't matter.
  X.S[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = embed, stdev = apply(embed, MARGIN = 2, FUN = sd),
    key = "PC_", assay = SeuratObject::DefaultAssay(X.S))
  X.S <- Seurat::FindNeighbors(
    X.S, reduction = "pca", dims = 1:n.dims,
    prune.SNN = prune.SNN, verbose = FALSE)

  if (sum(X.S$da.up) > 1) {
    up.S <- subset(X.S, subset = da.up == TRUE & da.down == FALSE)
    up.S <- Seurat::FindNeighbors(
      up.S, reduction = "pca", dims = 1:n.dims, verbose = FALSE)
    up.S <- Seurat::FindClusters(
      up.S, resolution = resolution,
      group.singletons = group.singletons,
      verbose = FALSE, ...)
    up.clusters <- as.numeric(up.S@active.ident)
    up.clusters[up.S@active.ident == "singleton"] <- 0
  } else {
    up.clusters <- NULL
  }

  n.up.clusters <- length(unique(up.clusters)) - as.numeric(0 %in% up.clusters)

  if (sum(X.S$da.down) > 1) {
    down.S <- subset(X.S, subset = da.down == TRUE & da.up == FALSE)
    down.S <- Seurat::FindNeighbors(
      down.S, reduction = "pca", dims = 1:n.dims, verbose = FALSE)
    down.S <- Seurat::FindClusters(
      down.S, resolution = resolution,
      group.singletons = group.singletons,
      verbose = FALSE, ...
    )
    down.clusters <- as.numeric(down.S@active.ident) + n.up.clusters
    down.clusters[down.S@active.ident == "singleton"] <- 0
  } else {
    down.clusters <- NULL
  }

  X.S$da.region.label <- 0
  X.S$da.region.label[X.S$da.up & !X.S$da.down] <- up.clusters
  X.S$da.region.label[X.S$da.down & !X.S$da.up] <- down.clusters

  lapply(obj_list, function(obj) {
    da_regions <- with(obj, {
      cell.names <- rownames(embeddings$pca)[cells$cell.idx]
      da.region.label <- X.S@meta.data[cell.names, "da.region.label"]
      cell.labels <- X.S@meta.data[cell.names, "sample"]

      significant <- logical(length(cell.names))
      significant[cells$da.up] <- TRUE
      significant[cells$da.down] <- TRUE
      da.region.label[!significant] <- 0

      # remove small clusters with cells < min.cell
      if (is.null(min.cell)) {
        min.cell <- as.integer(colnames(cells$da.ratio)[1])
        cat("Using min.cell = ", min.cell, "\n", sep = "")
      }

      da.region.label.tab <- table(da.region.label)
      if (min(da.region.label.tab) < min.cell) {
        da.region.to.remove <- as.numeric(names(da.region.label.tab)[
          which(da.region.label.tab < min.cell)])
        cat("Removing ", length(da.region.to.remove),
            " DA regions with cells < ", min.cell, ".\n", sep = "")
        da.region.label.old <- da.region.label
        for (ii in da.region.to.remove) {
          da.region.label[da.region.label.old == ii] <- 0
          # Do not change the name of DA cluster
          #da.region.label[da.region.label.old > ii] <- da.region.label[
          #  da.region.label.old > ii] - 1
        }
      }

      X.da <- unique(da.region.label)
      X.da <- X.da[X.da != 0]
      X.n.da <- length(X.da)
      X.da.stat <- matrix(0, nrow = X.n.da, ncol = 3)
      colnames(X.da.stat) <- c("DA.score", "pval.wilcoxon", "pval.ttest")
      rownames(X.da.stat) <- X.da
      if (X.n.da > 0) {
        for (ii in X.da) {
          X.da.stat[as.character(ii), ] <- DAseq:::getDAscore(
            cell.labels = cell.labels,
            cell.idx = which(da.region.label == ii),
            labels.1 = args$labels.1, labels.2 = args$labels.2
          )
        }
      } else {
        warning("No DA regions found.")
      }

      list(
        cell.idx = cells$cell.idx,
        da.region.label = da.region.label,
        DA.stat = X.da.stat,
        da.region.plot = NULL
      )
    })

    obj$regions <- da_regions
    obj
  })
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
#' for \code{Seurat} and \code{STG}.
#' @param return_raw Return raw results produced by specified \code{method}.
#' @param GPU Which GPU to use (GPU IDs), default using CPU. Note: this value
#' will be used to set \code{CUDA_VISIBLE_DEVICES} environment.
#' @param ... Additional arguments passed to marker finder.
#' @return A list of marker genes
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
      da <- unique(obj$da)
      da.regions.to.run <- da[da != 0]
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
      # FIXME: Do not use length to calculate DA subpopulations
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
      # FIXME: Do not use length to calculate DA subpopulations
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
#' @param reduction Cell embeddings to plot
#' @param order Order cells by DA score
#' @param label Label cell custers
#' @param circle Draw circle around clusters
#' @param point_size Set point size
#' @param point_alpha Set point alpha
#' @param text_size Set text size
#' @param theme_size Set theme size
#' @return ggplot object
#' @export
plotDACellScore <- function(
  obj, reduction,
  order = FALSE,
  label = FALSE,
  circle = FALSE,
  point_size = 1,
  point_alpha = 1,
  text_size = 4,
  theme_size = 10
) {
  embedding <- obj$embeddings[[reduction]][obj$cells$cell.idx, ]
  data_to_plot <- data.frame(
    Dim1 = embedding[, 1],
    Dim2 = embedding[, 2],
    clusters = obj$info$ident[obj$cells$cell.idx],
    score = obj$cells$da.pred
  )

  if (order)
    data_to_plot <- dplyr::arrange(data_to_plot, score)

  p <- data_to_plot %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(
      mapping = ggplot2::aes(x = Dim1, y = Dim2, col = score),
      size = point_size, alpha = point_alpha) +
    ggplot2::scale_color_gradientn(colours = c("blue", "white", "red")) +
    cowplot::theme_cowplot(theme_size) &
    ggplot2::theme(legend.title = ggplot2::element_blank())

  if (label) {
    label_positions <- data_to_plot %>%
      dplyr::group_by(clusters) %>%
      dplyr::summarise(x = median(Dim1), y = median(Dim2))

    p <- p +
      ggrepel::geom_text_repel(
        data = label_positions,
        mapping = ggplot2::aes(
          x = x, y = y,
          label = clusters
        ),
        size = text_size,
        show.legend = FALSE,
        bg.color = "white",
        fontface = "bold"
      )
  }

  if (circle) {
    p <- p + ggplot2::stat_ellipse(
      mapping = ggplot2::aes(x = Dim1, y = Dim2, group = clusters)
    )
  }

  p
}

#' Plot DA Embedding
#'
#' @inheritParams plotDACellScore
#' @export
plotDAEmbedding <- function(obj, reduction, group_by = "group") {
  embedding <- obj$embeddings[[reduction]]
  data_to_plot <- data.frame(
    Dim1 = embedding[, 1],
    Dim2 = embedding[, 2],
    Group = obj$info[[group_by]]
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
#' @param color_fun Color function to plot clusters.
#' Such as \code{\link[Seurat]{DiscretePalette}}
#' @export
plotDACellLabel <- function(obj, reduction, color_fun = scales::hue_pal()) {
  embedding <- obj$embeddings[[reduction]][obj$cells$cell.idx, ]
  X.da.label <- obj$regions$da.region.label
  X.da.order <- order(X.da.label, decreasing = FALSE)
  X.n.da <- length(unique(X.da.label)) - 1
  DAseq::plotCellLabel(
    X = embedding[X.da.order, ],
    label = as.factor(X.da.label[X.da.order])
  ) +
    ggplot2::scale_color_manual(
      values = c("gray", color_fun(X.n.da))
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
