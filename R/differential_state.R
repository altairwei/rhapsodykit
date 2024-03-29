#' Perform differential state analysis on pseudo-bulk data.
#'
#' @inheritParams calculate_pseudo_bulk
#' @param pb Object return by \code{\link{calculate_pseudo_bulk}}
#' @param contrasts Character vector specifying contrasts,
#'  see \code{\link[limma]{makeContrasts}}. If \code{NULL}, all possible
#'  combination will be tested.
#' @param method Specify which bulk RNA-seq DE methods to apply.
#' @param ... pass to \code{\link[muscat]{pbDS}}
#' @export
pseudobulk_diff_state <- function(
  sce, pb, contrasts = NULL,
  method = c("edgeR", "DESeq2", "limma-trend", "limma-voom"),
  ...
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
    ...
  )

  res
}

#' Find group-specific genes.
#'
#' @param target_group Which group to find marker genes.
#' @inheritParams pseudobulk_diff_state
#' @export
find_group_marker_genes <- function(
  sce, pb, target_group,
  method = c("edgeR", "DESeq2", "limma-trend", "limma-voom")
) {
  ei <- S4Vectors::metadata(sce)$experiment_info

  group_lvs <- levels(ei$group_id)
  group_lvs <- group_lvs[group_lvs != target_group]

  ei$group_id <- forcats::fct_collapse(
    ei$group_id, ..OTHERS = group_lvs)

  mm <- model.matrix(~ 0 + group_id, ei)
  dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))

  cm <- limma::makeContrasts(
    contrasts = paste(target_group, "..OTHERS", sep = "-"), levels = mm)

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
#' @inheritParams diff_state_filter
#' @param fdr_limit Upper limit of FDR.
#' @param logfc_limit Lower limit of absolute value of logFC.
#' @return Same as \code{tbl_list}
#' @export
diff_state_significant <- function(
  results, fdr_limit = 0.05, logfc_limit = 1) {
  results %>%
    diff_state_filter(
      p_adj.loc < fdr_limit,
      abs(logFC) > logfc_limit
    ) %>%
    diff_state_apply(dplyr::arrange, p_adj.loc)
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
diff_state_significant_lfc <- function(
  ds, contrasts, clusters, FDR = 0.05, logFC = 1) {
  df <- ds %>%
      diff_state_significant(FDR, logFC) %>%
      diff_state_pull(contrasts, clusters, c("gene", "logFC"))

  fc_list <- df[, "logFC"]
  names(fc_list) <- as.character(df[, "gene"])

  fc_list
}

#' Prepare ranked logFC vector named with gene names.
#'
#' @inheritParams diff_state_significant_lfc
#' @return A ranked logFC vector named gene names
#' @export
diff_state_ranked_lfc <- function(ds, contrasts, clusters) {
  df <-  ds %>%
    # We don't need to filter genes with given cutoff
    diff_state_pull(contrasts, clusters, c("gene", "logFC"))

  fc_list <- df[, "logFC"]
  names(fc_list) <- as.character(df[, "gene"])

  ranked_list <- sort(fc_list, decreasing = TRUE)

  ranked_list
}

#' Filter differential state analysis results.
#'
#' @inheritParams diff_state_apply
#' @param ... Arguments will be passed to \code{\link[dplyr]{filter}}
#' @return Same as \code{tbl_list}
#' @export
diff_state_filter <- function(results, ...) {
  diff_state_apply(results, dplyr::filter, ...)
}

#' Apply funtion on each data.frame of differential state analysis results
#' @param results Results returned by \code{\link[muscat]{pbDS}} or
#'  \code{\link{pseudobulk_diff_state}}
#' @param ... Arguments will be passed to \code{fun}
#' @export
diff_state_apply <- function(results, fun, ...) {
  stopifnot(check_diff_state_results(results))
  results$table <- lapply(results$table, function(contrast) {
    lapply(contrast, function(subpopulation) fun(subpopulation, ...))
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


#' Pull out column from diff state results.
#'
#' @inheritParams diff_state_apply
#' @param contrasts choose a contrast.
#' @param clusters choose a cluster.
#' @param columns choose a column.
#' @export
diff_state_pull <- function(
  results, contrasts, clusters, columns, simplify = TRUE) {
  stopifnot(check_diff_state_results(results))
  res <- lapply(results$table[contrasts], function(contrast) {
    lapply(contrast[clusters], function(subpopulation) subpopulation[columns])
  })

  if (simplify) {
    if (length(contrasts) == 1) {
      if (length(clusters) == 1) {
        if (length(columns) == 1) {
          # a vector
          return(res[[1]][[1]][[1]])
        } else {
          # a data.frame
          return(res[[1]][[1]])
        }
      } else {
        if (length(columns) == 1) {
          # a vector under each cluster
          lapply(res[[1]], function(x) x[[1]])
        } else {
          # a data.frame under each cluster
          return(res[[1]])
        }
      }
    } else {
      if (length(clusters) == 1) {
        if (length(columns) == 1) {
          # return a vector under each contrast
          lapply(res, function(x) x[[1]][[1]])
        } else {
          # return a data.frame on each contrast
          lapply(res, function(x) x[[1]])
        }
      } else {
        if (length(columns) == 1) {
          # vector under 2-level nested list
          lapply(res, function(x) {
            lapply(x, function(y) y[[1]])
          })
        } else {
          # data.frame under 2-level nested list
          return(res)
        }
      }
    }
  } else {
    return(res)
  }
}

#' Subset Differential State Results
#'
#' @inheritParams diff_state_pull
#' @export
diff_state_subset <- function(results, contrasts = NULL, clusters = NULL) {
  if (!is.null(contrasts))
    results$table <- results$table[contrasts]

  if (!is.null(clusters)) {
    results$table <- purrr::map(results$table, function(contr) contr[clusters])
    results$data <- results$data[clusters]
    results$fit <- results$fit[clusters]
  }

  results
}