#' Perform differential state analysis on pseudo-bulk data.
#'
#' @inheritParams calculate_pseudo_bulk
#' @param pb Object return by \code{\link{calculate_pseudo_bulk}}
#' @param contrasts Character vector specifying contrasts,
#'  see \code{\link[limma]{makeContrasts}}. If \code{NULL}, all possible
#'  combination will be tested.
#' @param method Specify which bulk RNA-seq DE methods to apply.
#' @export
pseudobulk_diff_state <- function(
  sce, pb, contrasts = NULL,
  method = c("edgeR", "DESeq2", "limma-trend", "limma-voom")
) {
  ei <- S4Vectors::metadata(sce)$experiment_info

  mm <- model.matrix(~ 0 + group_id, ei)
  dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))

  if (is.null(contrasts)) {
    contrasts <- sapply(
      combn(rev(levels(ei$group_id)), 2, simplify = FALSE),
      function(x) paste(x, collapse = "-"))
  }

  cm <- limma::makeContrasts(contrasts = contrasts, levels = mm)

  res <- muscat::pbDS(
    pb,
    method = match.arg(method),
    design = mm,
    contrast = cm,
    verbose = FALSE
  )

  res$table
}

#' Filter differential state analysis results by FDR and logFC.
#'
#' @param tbl_list A list of data.frame that are one of results from
#'  \code{\link{pseudobulk_diff_state}}
#' @param fdr_limit Upper limit of FDR.
#' @param logfc_limit Lower limit of absolute value of logFC.
#' @return Same as \code{tbl_list}
#' @export
diff_state_filter <- function(
  tbl_list,
  fdr_limit = 0.05,
  logfc_limit = 1
) {
  stopifnot(all(sapply(tbl_list, is.data.frame)))

  tbl_filtered <- lapply(tbl_list, function(subpopulation) {
    subpopulation <- dplyr::filter(
      subpopulation,
      p_adj.loc < fdr_limit,
      abs(logFC) > logfc_limit
    )
    dplyr::arrange(subpopulation, p_adj.loc)
  })

  tbl_filtered
}

#' Filter differential state analysis results by gene expression frequencies.
#'
#' Only retain genes that are expressed in an average of given percentage of
#'  cells in at least 1 group.
#'
#' @inheritParams diff_state_filter
#' @inheritParams calculate_pseudo_bulk
#' @param percent Minimum percentage of cells in which a gene expressed.
#' @inherit diff_state_filter return
#' @export
expr_freq_filter <- function(tbl_list, sce, percent = 0.1) {
  # one assay of dimensions #genes x #samples/#groups per cluster.
  frq <- muscat::calcExprFreqs(sce, assay = "counts")

  # Get genes that are expressed in an average of `percent` of cells in
  # at least 1 group.
  group_ids <- levels(sce$group_id)
  frq10 <- vapply(
    as.list(SummarizedExperiment::assays(frq)),
    function(cluster_frq) apply(cluster_frq[, group_ids] > percent, 1, any),
    logical(nrow(sce))
  )

  cluster_ids <- levels(sce$cluster_id)
  names(cluster_ids) <- cluster_ids

  # Retain genes that are passed filter criteria
  tbl_list2 <- lapply(cluster_ids, function(k) {
    dplyr::filter(
      tbl_list[[k]],
      gene %in% names(which(frq10[, k]))
    )
  })

  tbl_list2
}

#' Summarize results from differential state analysis.
#'
#' @inheritParams diff_state_filter
#' @inheritParams calculate_pseudo_bulk
#' @export
diff_state_summary <- function(tbl_list, sce) {
  # nb. of DS genes & % of total by cluster
  n_de <- vapply(tbl_list, nrow, numeric(1))
  p_de <- format(n_de / nrow(sce) * 100, digits = 3)

  data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)
}

#' Format differential state testing results into an easily filterable table
#'
#' @inheritParams diff_state_filter
#' @export
diff_state_format <- function(tbl_list) {
  dplyr::bind_rows(tbl_list)
}
