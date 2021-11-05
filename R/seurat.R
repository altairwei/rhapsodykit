#' Generate plots with the Seurat object
#'
#' @param object A seurat object.size.
#' @param output_folder Folder to place plots.
#'
#' @export
generate_seurat_plots <- function(object, output_folder) {
  # Quality Control
  tryCatch(
    expr = {
      p_qc_metrics <- Seurat::VlnPlot(
        object,
        features = c("nFeature_RNA", "nCount_RNA"),
        ncol = 2
      )
      p_qc_metrics_image_width <- 6

      if (!is.null(object@meta.data[["percent.mt"]])) {
        p_qc_metrics <- p_qc_metrics + Seurat::VlnPlot(
          object, features = "percent.mt")
        p_qc_metrics_image_width <- p_qc_metrics_image_width + 3
      }

      if (!is.null(object@meta.data[["percent.cp"]])) {
        p_qc_metrics <- p_qc_metrics + Seurat::VlnPlot(
          object, features = "percent.cp")
        p_qc_metrics_image_width <- p_qc_metrics_image_width + 3
      }
      save_plot(
        file.path(output_folder, "QC_Metrics.png"),
        p_qc_metrics, width = p_qc_metrics_image_width
      )
      save_plot(
        file.path(output_folder, "QC_Count_Feature_Scatter.png"),
        Seurat::FeatureScatter(
          object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
      )
    },
    error = function(e) message(toString(e))
  )

  # Feature selection
  tryCatch(
    expr = {
      # Identify the 10 most highly variable genes
      feature_selection_top10 <- head(
        Seurat::VariableFeatures(object), 10)
      # plot variable features with and without labels
      p_feature_selection_1 <- Seurat::VariableFeaturePlot(object)
      p_feature_selection_2 <- Seurat::LabelPoints(
        plot = p_feature_selection_1,
        points = feature_selection_top10, repel = TRUE)

      save_plot(
        file.path(output_folder, "Feature_Selection.png"),
        p_feature_selection_1 + p_feature_selection_2,
        width = 12
      )
    },
    error = function(e) message(toString(e))
  )

  # PCA Plots
  tryCatch(
    expr = {
      npc <- length(object[["pca"]]@stdev)
      p_pca_dim_load <- Seurat::VizDimLoadings(
        object, dims = 1:4, reduction = "pca")
      save_plot(
        file.path(output_folder, "PCA_Dim_Loadings.png"),
        p_pca_dim_load
      )
      p_pca_dim_scatter <- Seurat::DimPlot(object, reduction = "pca")
      save_plot(
        file.path(output_folder, "PCA_Scatter.png"),
        p_pca_dim_scatter
      )
    },
    error = function(e) message(toString(e))
  )

  tryCatch(
    expr = {
      npc <- length(object[["pca"]]@stdev)
      p_jackstraw <- Seurat::JackStrawPlot(object, dims = 1:npc) +
        ggplot2::theme(legend.position = "none")
      p_elbow <- Seurat::ElbowPlot(object, ndims = npc)
      save_plot(
        file.path(output_folder, "PCA_Dimensionality.png"),
        p_jackstraw + p_elbow,
        width = 12
      )
    },
    error = function(e) message(toString(e))
  )

  tryCatch(
    expr = {
      p_pca_dim_heatmap <- Seurat::DimHeatmap(
        object, dims = 1:9, cells = 500, balanced = TRUE)
      save_plot(
        file.path(options$output_folder, "PCA_Dim_Heatmap.png"),
        p_pca_dim_heatmap,
        width = 12
      )
    },
    error = function(e) message(toString(e))
  )

  # UMAP/tSNE Plots
  if (Seurat::DefaultAssay(object) == "integrated") {
    tryCatch(
      expr = {
        p_umap_sample <- Seurat::DimPlot(
          object, reduction = "umap", group.by = "sample")
        p_umap_cluster <- Seurat::DimPlot(
          object, reduction = "umap", label = TRUE)
        save_plot(
          file.path(output_folder, "UMAP_Scatter.png"),
          p_umap_sample + p_umap_cluster,
          width = 12
        )

        p_tsne_sample <- Seurat::DimPlot(
          object, reduction = "tsne", group.by = "sample")
        p_tsne_cluster <- Seurat::DimPlot(
          object, reduction = "tsne", label = TRUE)
        save_plot(
          file.path(output_folder, "TSNE_Scatter.png"),
          p_tsne_sample + p_tsne_cluster,
          width = 12
        )


        sample <- unique(unlist(object[["sample"]]))
        n_sample <-  length(sample)
        p_umap_split_sample <- Seurat::DimPlot(
          object, reduction = "umap", split.by = "sample")
        p_tsne_split_sample <- Seurat::DimPlot(
          object, reduction = "tsne", split.by = "sample")
        save_plot(
          file.path(output_folder, "UMAP_Scatter_By_Sample.png"),
          p_umap_split_sample, width = 6 * n_sample
        )
        save_plot(
          file.path(output_folder, "TSNE_Scatter_By_Sample.png"),
          p_tsne_split_sample, width = 6 * n_sample
        )
      },
      error = function(e) message(toString(e))
    )
  } else {
    tryCatch(
      expr = {
        p_umap_dim <- Seurat::DimPlot(
          object, reduction = "umap", label = TRUE)
        save_plot(
          file.path(output_folder, "UMAP_Scatter.png"),
          p_umap_dim
        )
        p_tsne_dim <- Seurat::DimPlot(
          object, reduction = "tsne", label = TRUE)
        save_plot(
          file.path(output_folder, "TSNE_Scatter.png"),
          p_tsne_dim
        )
      },
      error = function(e) message(toString(e))
    )
  }

}

#' Find gene markers for all cell identity classes
#'
#' @param object A Seurat object.
#' @param output_folder Folder to place outputs.
#'
#' @export
perform_find_all_markers <- function(object, output_folder) {
  markers_df <- Seurat::FindAllMarkers(
    object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  readr::write_csv(markers_df, file.path(output_folder, "Markers_All.csv"))
  top10 <- dplyr::top_n(
    dplyr::group_by(markers_df, cluster), n = 10, wt = avg_log2FC)
  marker_heatmap <- Seurat::DoHeatmap(
    object, features = top10$gene) + Seurat::NoLegend()
  save_plot(
    file.path(output_folder, "Markers_Top10_Heatmap.png"),
    marker_heatmap, width = 14, height = 14)
}

#' Perform single sample Seurat analysis
#'
#' @param object A Seurat object
#' @param mt_gene_file A file contains mitochondrial genes.
#' @param cp_gene_file A file contains chloroplast genes.
#' @param dimensionality Number of dimensions to use as input.
#' @return A Seurat object
#'
#' @export
single_sample_analysis <- function(
  object,
  mt_gene_file = NULL,
  cp_gene_file = NULL,
  dimensionality = 20
) {
  # Quality Control
  if (!is.null(mt_gene_file)) {
    mt_gene_list <- readr::read_lines(mt_gene_file)
    mt_gene_list <- intersect(rownames(object[["RNA"]]), mt_gene_list)
    if (length(mt_gene_list) > 0) {
      message(sprintf(
        "Found %d mitochondrial genes in expression matrix.",
        length(mt_gene_list)))
      object[["percent.mt"]] <- Seurat::PercentageFeatureSet(
        object, features = mt_gene_list)
    } else {
      message("No mitochondrial genes were found.")
    }
  }

  if (!is.null(cp_gene_file)) {
    cp_gene_list <- readr::read_lines(cp_gene_file)
    cp_gene_list <- intersect(rownames(object[["RNA"]]), cp_gene_list)
    if (length(cp_gene_list) > 0) {
      message(sprintf(
        "Found %d chloroplast genes in expression matrix.",
        length(mt_gene_list)))
      object[["percent.cp"]] <- Seurat::PercentageFeatureSet(
        object, features = cp_gene_list)
    } else {
      message("No chloroplast genes were found.")
    }
  }

  # TODO: 增加 SCTransform 的选项！
  # Normalizing the data
  seurat_obj <- Seurat::NormalizeData(
    object, normalization.method = "LogNormalize", scale.factor = 10000)

  # Feature selection
  seurat_obj <- Seurat::FindVariableFeatures(
    seurat_obj, selection.method = "vst", nfeatures = 2000)

  # Scaling the data
  all_genes <- rownames(seurat_obj)
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = all_genes)

  # PCA Plots

  # Perform linear dimensional reduction
  seurat_obj <- Seurat::RunPCA(
    seurat_obj, features = Seurat::VariableFeatures(seurat_obj))
  npc <- length(seurat_obj[["pca"]]@stdev)

  # Determine the 'dimensionality' of the dataset
  seurat_obj <- Seurat::JackStraw(seurat_obj, dims = npc, num.replicate = 100)
  seurat_obj <- Seurat::ScoreJackStraw(seurat_obj, dims = 1:npc)

  # Cluster the cells
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:dimensionality)
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.5)

  # UMAP/tSNE Plots
  #TODO: Make sure umap-learn work properly
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:dimensionality)
  #TODO: Check duplicates manually
  # Workaround: https://github.com/satijalab/seurat/issues/167
  seurat_obj <- Seurat::RunTSNE(
    seurat_obj, dims = 1:dimensionality, check_duplicates = FALSE)

  seurat_obj
}

#' Merge Seurat objects into one without batch effect correction
#'
#' @inheritParams integrated_sample_analysis
#' @param normalization Normalization method
#' @return A merged Seurat object
#'
#' @export
merge_sample <- function(
  obj_list, normalization = "LogNormalize", analyze = TRUE
) {
  stopifnot(all(sapply(obj_list, inherits, "Seurat")))
  stopifnot(length(obj_list) >= 2)

  cell_id_prefix <- NULL
  if (!is.null(names(obj_list))) {
    cell_id_prefix <- names(obj_list)
  }

  obj_merged <- merge(
    x = obj_list[[1]],
    y = obj_list[-1],
    add.cell.ids = cell_id_prefix
  )

  if (analyze)
    obj_merged <- merge_analysis(obj_merged, normalization = normalization)

  obj_merged
}

merge_analysis <- function(
  object, n_dims = 20, cluster_res = 0.5, normalization = "LogNormalize"
) {

  if (normalization == "SCTransform") {
    object <- Seurat::SCTransform(
      object, do.scale = FALSE)
  } else {
    object <- Seurat::NormalizeData(
      object, normalization.method = normalization)
  }

  object <- Seurat::FindVariableFeatures(
      object, selection.method = "vst", nfeatures = 2000)

  # Run the standard workflow for visualization and clustering
  object <- Seurat::ScaleData(object)
  object <- Seurat::RunPCA(object)

  object <- Seurat::RunUMAP(object, dims = 1:n_dims)

  #TODO: Check duplicates manually
  # Workaround: https://github.com/satijalab/seurat/issues/167
  object <- Seurat::RunTSNE(
    object, dims = 1:n_dims, check_duplicates = FALSE)

  object <- Seurat::FindNeighbors(object, dims = 1:n_dims)
  object <- Seurat::FindClusters(object, resolution = cluster_res)

  object
}

#' Perform Seurat integrated analysis for grouped samples
#'
#' @inheritParams integration_anchorset
#' @param analyze Perform general Seurat analysis on integrated object or not.
#' @param ... pass to \code{integration_anchorset}
#' @return A integrated Seurat object.
#'
#' @export
integrated_sample_analysis <- function(
  obj_list, n_dims = 20, analyze = TRUE, ...
) {
  obj_anchors <- integration_anchorset(
    obj_list = obj_list, n_dims = n_dims, ...)

  obj_combined <- Seurat::IntegrateData(
    anchorset = obj_anchors, dims = 1:n_dims)

  if (analyze)
    obj_combined <- integration_analysis(obj_combined)

  obj_combined
}

#' Get integration anchorset from Seurat objects.
#'
#' @param obj_list A list of Seurat objects.
#' @param n_dims Number of dimensions to use as input.
#' @param normalization Normalization method.
#' @param ... pass to \code{\link[Seurat]{FindIntegrationAnchors}
#' @inheritParams Seurat::FindIntegrationAnchors
#' @return Returns an \code{AnchorSet} object that can be used as input
#'  to \code{\link[Seurat]{IntegrateData}}.
#'
#' @export
integration_anchorset <- function(
  obj_list,
  n_dims = 20,
  reduction = c("cca", "rpca"),
  k.anchor = 5,
  reference = NULL,
  normalization = "LogNormalize",
  ...
) {
  stopifnot(all(sapply(obj_list, inherits, "Seurat")))

  obj_list <- lapply(obj_list, function(obj) {

    if (normalization == "SCTransform") {
      obj <- Seurat::SCTransform(obj, do.scale = FALSE)
    } else {
      obj <- Seurat::NormalizeData(obj, normalization.method = normalization)
    }

    obj <- Seurat::FindVariableFeatures(
      obj, selection.method = "vst", nfeatures = 2000)
  })

  features <- Seurat::SelectIntegrationFeatures(object.list = obj_list)

  reduction <- match.arg(reduction)
  if (reduction == "rpca") {
    obj_list <- lapply(obj_list, function(x) {
        x <- Seurat::ScaleData(x, features = features, verbose = FALSE)
        x <- Seurat::RunPCA(x, features = features, verbose = FALSE)
    })
  }

  obj_anchors <- Seurat::FindIntegrationAnchors(
    obj_list,
    dims = 1:n_dims,
    reduction = reduction,
    anchor.features = features,
    k.anchor = k.anchor,
    reference = reference,
    ...
  )

  obj_anchors
}

#' Perform general Seurat analysis on integrated object
#'
#' @param object Integrated Seurat object
#' @param n_dims Number of dimensions to use as input.
#' @param cluster_res Resolution for \code{\link[Seurat]{FindClusters}}
#' @return A Seurat object
#' @export
integration_analysis <- function(
  object, n_dims = 20, cluster_res = 0.5
) {
  Seurat::DefaultAssay(object) <- "integrated"

  # Run the standard workflow for visualization and clustering
  object <- Seurat::ScaleData(object)
  object <- Seurat::RunPCA(object)
  # t-SNE and Clustering

  #TODO: Make sure umap-learn work properly
  object <- Seurat::RunUMAP(object, dims = 1:n_dims)

  #TODO: Check duplicates manually
  # Workaround: https://github.com/satijalab/seurat/issues/167
  object <- Seurat::RunTSNE(
    object, dims = 1:n_dims, check_duplicates = FALSE)

  object <- Seurat::FindNeighbors(object, dims = 1:n_dims)
  object <- Seurat::FindClusters(object, resolution = cluster_res)

  object
}

#' Find conserved markers for all cell identity classes
#'
#' @param object A Seurat object.
#' @return A data.frame containing conserved markers.
#'
#' @export
find_all_conserved_markers <- function(object) {
  Seurat::DefaultAssay(object) <- "RNA"
  idents_all <- sort(unique(Seurat::Idents(object)))
  df_list <- lapply(idents_all, function(i) {
    df <- tryCatch(
      Seurat::FindConservedMarkers(
        object, ident.1 = i, grouping.var = "sample"),
      error = function(e) {
        message(toString(e))
        NULL
      }
    )
    if (!is.null(df)) {
      df[["cluster"]] <- i
      df[["gene"]] <- rownames(df)
      rownames(df) <- NULL
      return(df)
    } else {
      return(NULL)
    }
  })
  do.call(dplyr::bind_rows, df_list)
}

#' Generate average expression of different sample
#'  for all cell identity classes
#'
#' @param object A Seurat object.
#' @return A data.frame containing average expression.
#'
#' @export
find_all_avg_expr_genes <- function(object) {
  idents_all <- sort(unique(Seurat::Idents(object)))
  df_list <- lapply(idents_all, function(i) {
    ident_cells <- subset(object, idents = i)
    # Specify identity of cells based on value of meta.data[["sample"]]
    Seurat::Idents(ident_cells) <- "sample"
    # AverageExpression will be applied to every `Assay` object.
    avg_ident_cells <- as.data.frame(
      log1p(Seurat::AverageExpression(ident_cells)$RNA))
    avg_ident_cells[["cluster"]] <- i
    avg_ident_cells[["gene"]] <- rownames(avg_ident_cells)
    avg_ident_cells
  })
  do.call(dplyr::bind_rows, df_list)
}

#' Find differential expressed genes on variable conditions
#'  for each cell cluster.
#'
#' @param object Seurat object.
#' @return A data frame.
#'
#' @export
find_all_diff_expr_genes <- function(object) {
  idents_all <- sort(unique(Seurat::Idents(object)))
  sample <- unique(unlist(object[["sample"]]))
  # Produce an permutation of all conditions
  all_comb <- as.data.frame(combn(sample, 2), stringsAsFactors = FALSE)
  names(all_comb) <- NULL
  # Find diff genes on each comparision
  comp_list <- lapply(all_comb, function(couple) {
    message(sprintf("Calculating comparision %s vs %s", couple[2], couple[1]))
    # Perform on each cell cluster
    diff_list <- lapply(idents_all, function(id) {
      ident_cells <- subset(object, idents = id)
      # Specify identity of cells based on value of meta.data[["sample"]]
      Seurat::Idents(ident_cells) <- "sample"
      message("Calculating cluster ", id)
      # We assume that the control group is the first one
      #   and the stimulated group is the second one.
      tryCatch(
        {
          diff_genes <- Seurat::FindMarkers(
            ident_cells, ident.1 = couple[2], ident.2 = couple[1])
          diff_genes[["test_group"]] <- couple[2]
          diff_genes[["ctrl_group"]] <- couple[1]
          diff_genes[["comparison"]] <- sprintf(
            "%s_vs_%s", couple[2], couple[1])
          diff_genes[["cluster"]] <- id
          diff_genes[["gene"]] <- rownames(diff_genes)
          diff_genes
        },
        error = function(e) {
          message(toString(e))
          NULL
        }
      )
    })
    do.call(dplyr::bind_rows, diff_list)
  })
  do.call(dplyr::bind_rows, comp_list)
}

#' Plot scatter diagram to compare average expression value of genes on
#'   variable conditions for each cell cluster.
#'
#' @param avg_genes A data frame produced by \code{find_all_avg_expr_genes}
#' @param diff_genes A data frame produced by \code{find_all_diff_expr_genes}
#' @return A named list of ggplot2 object.
#'
#' @export
plot_avg_expr_genes <- function(avg_genes, diff_genes) {
  diff_genes_by_comp <- split(diff_genes, diff_genes[["comparison"]])
  results <- lapply(diff_genes_by_comp, function(diff_genes_comp) {
    # test condition should be placed on y-axis, ctrl condition on x-axis
    condition_ctrl <- unique(diff_genes_comp[["ctrl_group"]])
    condition_test <- unique(diff_genes_comp[["test_group"]])
    # Calculate coefficient of correlation
    #TODO: fill NA with zero.
    df_by_cluster <- split(avg_genes, avg_genes[["cluster"]])
    xmin <- min(avg_genes[[condition_ctrl]], na.rm = TRUE)
    ymax <- max(avg_genes[[condition_test]], na.rm = TRUE)
    df_by_cluster <- lapply(df_by_cluster, function(x) {
      label_df <- tryCatch(
        {
          pearson_r <- cor(x[[condition_ctrl]], x[[condition_test]],
            method = "pearson", use = "complete.obs")
          data.frame(
            coeff_label = sprintf(
              "'Pearson ' * italic(R^2) == '%.2f'", pearson_r^2),
            cluster = unique(x[["cluster"]]),
            position_x = xmin,
            position_y = ymax
          )
        },
        error = function(e) {
          message(
            "Error occurred when processing cluster ",
            unique(x[["cluster"]]), " of ",
            condition_test, "_vs_", condition_ctrl)
          message(toString(e))
          data.frame(
            coeff_label = "",
            cluster = unique(x[["cluster"]]),
            position_x = xmin,
            position_y = ymax
          )
        }
      )
      label_df
    })
    df_coff_label <- do.call(dplyr::bind_rows, df_by_cluster)
    # Highlight DEGs
    diff_genes_to_highlight <- dplyr::filter(diff_genes_comp, p_val < 0.0001)
    df_to_plot <- avg_genes
    df_to_plot$highlight <- "no"
    df_to_plot$highlight[
      df_to_plot$gene %in% diff_genes_to_highlight$gene] <- "yes"
    # Plot facets
    p1 <- ggplot2::ggplot(df_to_plot, ggplot2::aes(
        x = !!ggplot2::sym(condition_ctrl), y = !!ggplot2::sym(condition_test),
        color = highlight, alpha = highlight)) +
      ggplot2::ggtitle(unique(diff_genes_comp[["comparison"]])) +
      ggplot2::geom_point() +
      ggplot2::geom_text(
        data = df_coff_label,
        # Don't inherit color or alpha
        inherit.aes = FALSE,
        mapping = ggplot2::aes(
          x = position_x, y = position_y, label = coeff_label),
        hjust = 0, vjust = 1, size = 8, parse = TRUE) +
      ggplot2::scale_color_manual(values = c("black", "red"), guide = FALSE) +
      ggplot2::scale_alpha_manual(values = c(0.05, 1), guide = FALSE) +
      ggplot2::facet_wrap(ggplot2::vars(cluster)) +
      cowplot::panel_border() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    # Plot whole correlation diagram
    pearson_r_whole <- cor(
      avg_genes[[condition_ctrl]], avg_genes[[condition_test]],
      method = "pearson", use = "complete.obs")
    label_whole <- data.frame(
      coeff_label = sprintf(
        "'Pearson ' * italic(R^2) == '%.2f'", pearson_r_whole^2),
      position_x = xmin,
      position_y = ymax
    )
    p2 <- ggplot2::ggplot(avg_genes, ggplot2::aes(
        x = !!ggplot2::sym(condition_ctrl),
        y = !!ggplot2::sym(condition_test))) +
      ggplot2::ggtitle(unique(diff_genes_comp[["comparison"]])) +
      ggplot2::geom_point() +
      ggplot2::geom_text(
        data = label_whole,
        # Don't inherit color or alpha
        inherit.aes = FALSE,
        mapping = ggplot2::aes(
          x = position_x, y = position_y, label = coeff_label),
        hjust = 0, vjust = 1, size = 8, parse = TRUE) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    patchwork::wrap_plots(p1, p2)
  })

  results
}

#' Perform differentially expressed genes analysis.
#'
#' @param object A Seurat object
#' @param output_folder Folder to place outputs
#' @param draw_plot Draw plots or not.
#'
#' @export
perform_diff_gene <- function(object, output_folder, draw_plot = TRUE) {
  tryCatch(
    expr = {
      avg_genes_df <- find_all_avg_expr_genes(object)
      readr::write_csv(
        avg_genes_df,
        file.path(output_folder, "DEG_Average_Expression.csv"))
      diff_genes_df <- find_all_diff_expr_genes(object)
      readr::write_csv(
        diff_genes_df, file.path(output_folder, "DEG_All.csv"))
      if (isTRUE(draw_plot)) {
        all_p_list <- plot_avg_expr_genes(avg_genes_df, diff_genes_df)
        for (p_name in names(all_p_list)) {
          save_plot(
            file.path(output_folder,
              sprintf("DEG_Avg_Expr_Scatter_%s.png", p_name)),
            all_p_list[[p_name]],
            width = 32, height = 16
          )
        }
      }
    },
    error = function(e) message(toString(e))
  )
}

#' Reassign sample name.
#'
#' @param object A Seurat object
#' @param sample_map A named character vecter, such as `c("old" = "new")`.
#'
#' @export
reassign_sample <- function(object, sample_map) {
  if (is.null(names(sample_map)) || !is.character(sample_map))
    stop("sample_map must be a named character vecter.")

  object$sample <- key_mapping(object$sample, sample_map)

  object
}