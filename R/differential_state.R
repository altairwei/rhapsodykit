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

  res
}

check_diff_state_results <- function(x) {
  all(c("table", "data", "fit", "args") %in% names(x))
}

#' Filter differential state analysis results by FDR and logFC.
#'
#' @param results Results returned by \code{\link[muscat]{pbDS}} or
#'  \code{\link{pseudobulk_diff_state}}
#' @param fdr_limit Upper limit of FDR.
#' @param logfc_limit Lower limit of absolute value of logFC.
#' @return Same as \code{tbl_list}
#' @export
diff_state_filter <- function(
  results,
  fdr_limit = 0.05,
  logfc_limit = 1
) {
  stopifnot(check_diff_state_results(results))

  results$table <- lapply(results$table, function(contrast) {
    lapply(contrast, function(subpopulation) {
      subpopulation <- dplyr::filter(
        subpopulation,
        p_adj.loc < fdr_limit,
        abs(logFC) > logfc_limit
      )
      dplyr::arrange(subpopulation, p_adj.loc)
    })
  })

  results
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
expr_freq_filter <- function(results, sce, percent = 0.1) {
  stopifnot(check_diff_state_results(results))

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
  results$table <- lapply(results$table, function(contrast) {
    lapply(cluster_ids, function(k) {
      dplyr::filter(
        contrast[[k]],
        gene %in% names(which(frq10[, k]))
      )
    })
  })

  results
}

#' Summarize results from differential state analysis.
#'
#' @inheritParams diff_state_filter
#' @export
diff_state_summary <- function(results) {
  stopifnot(check_diff_state_results(results))

  purrr::imap(results$table, function(contrast, key) {
    # nb. of DS genes & % of total by cluster
    n_de <- vapply(contrast, nrow, numeric(1), USE.NAMES = FALSE)
    p_de <- format(n_de / nrow(results$args$pb) * 100, digits = 3)

    data.frame(
      "cluster_id" = names(contrast),
      "num_of_DS" = n_de,
      "percent_of_DS" = p_de,
      "contrast" = key
    )
  }) %>% dplyr::bind_rows()
}

#' Format differential state testing results into an easily filterable table
#'
#' @inheritParams diff_state_filter
#' @export
diff_state_format <- function(results) {
  results$table %>%
    lapply(dplyr::bind_rows) %>%
    dplyr::bind_rows()
}