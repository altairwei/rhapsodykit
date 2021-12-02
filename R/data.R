#' Wrapper for \code{ggsave} with more messages.
#'
#' @inheritParams ggplot2::ggsave
#'
#' @export
save_plot <- function(filename, plot, width = 7, height = 7, ...) {
  message(sprintf("Saving %d x %d in image: %s", width, height, filename))
  ggplot2::ggsave(filename, plot, width = width, height = height, ...)
}

save_png <- function(filename, plot, width = 7, height = 7, dpi = 300) {
  message(sprintf("Saving %d x %d in image: %s", width, height, filename))
  png(filename, width = width, height = height, units = "cm", res = dpi)
  print(plot)
  dev.off()
}

#' Read raw expression matrix from Rhapsody WTA pipeline
#'
#' @param base_dir Root folder of Rhapsody WTA pipeline results
#' @return An object constructed from \code{\link[Matrix]{Matrix}}
#' @export
read_raw_csv <- function(base_dir = ".") {
  if (!dir.exists(base_dir)) {
    stop("Directory provided does not exist")
  }

  matrix_loc <- Sys.glob(file.path(base_dir, "*_RSEC_MolsPerCell.csv"))

  if (length(matrix_loc) != 1) {
    stop("`*_RSEC_MolsPerCell.csv` missing or more than one file was found.")
  }

  message(sprintf("Reading %s", matrix_loc))
  df <- readr::read_csv(matrix_loc, comment = "#", progress = TRUE,
    col_types = readr::cols(
      Cell_Index = readr::col_character()
    )
  )

  message("Matrix constructing")
  m <- data.matrix(df[, -1])
  rownames(m) <- df[["Cell_Index"]]
  m <- t(m)
  mtx <- Matrix::Matrix(m)

  rm(df, m)
  gc()

  mtx
}

#' Save Matrix data to disk.
#'
#' @param mtx An \code{\link[Matrix]{Matrix}} object produced
#'  from \code{read_raw_csv}
#' @param mtx_file Output filename
#' @export
save_expression_matrix <- function(mtx, mtx_file) {
  colnames_file <- sprintf("%s.colnames", mtx_file)
  rownames_file <- sprintf("%s.rownames", mtx_file)

  message(sprintf("Writting %s", mtx_file))
  Matrix::writeMM(mtx, file = mtx_file)
  write(colnames(mtx), colnames_file)
  write(rownames(mtx), rownames_file)
}

#' Read \code{Matrix} data from cached folder.
#'
#' @param base_dir Folder contains cached \code{\link[Matrix]{Matrix}} object.
#' @return An object constructed from \code{\link[Matrix]{Matrix}}
#' @export
read_mtx <- function(base_dir = ".") {
  if (!dir.exists(base_dir)) {
    stop("Directory provided does not exist")
  }

  matrix_loc <- Sys.glob(
    file.path(base_dir, "*_Expression_Matrix.mtx"))

  if (length(matrix_loc) != 1) {
    stop("`*_Expression_Matrix.mtx` missing or more than one file was found.")
  }

  colnames_loc <- sprintf("%s.colnames", matrix_loc)
  rownames_loc <- sprintf("%s.rownames", matrix_loc)
  stopifnot(file.exists(colnames_loc, rownames_loc))

  message(sprintf("Reading %s", matrix_loc))
  m <- Matrix::readMM(matrix_loc)
  colnames(m) <- readr::read_lines(colnames_loc)
  rownames(m) <- readr::read_lines(rownames_loc)

  m
}

#' Read expression data from results folder of Rhapsody WTA pipeline.
#'
#' @param base_dir The path of results folder.
#' @param use_mtx Specify whether to use the sparse matrix constructed
#'    by the \code{Matrix} package.
#' @return A matrix whose column names represent cell index,
#'    row names represent genes.
#'
#' @export
read_rhapsody_wta <- function(base_dir, use_mtx = FALSE) {
  expr_matrix <- NULL

  if (isTRUE(use_mtx)) {
    expr_matrix <- read_mtx(base_dir)
  } else {
    expr_matrix <- read_raw_csv(base_dir)
  }

  expr_matrix
}

#' Save Single-Cell Object to Disk
#'
#' This function can be used to save single-cell data object to
#' disk file.
#'
#' @param obj single-cell data object
#' @param filename output filename.
#' @param ... additional parameters
#' @export
save_to_disk <- function(obj, filename, ...) {
  UseMethod("save_to_disk")
}

#' @describeIn save_to_disk Save Seurat Object to Disk
#' @method save_to_disk Seurat
#' @export
save_to_disk.Seurat <- function(obj, filename, ...) {
  SeuratDisk::SaveH5Seurat(obj, filename, overwrite = TRUE)
}

#' Access Cellular Data from SeuratDisk File
#'
#' Retrieves data (feature expression, PCA scores, metrics, etc.) for a set
#' of cells in a SeuratDisk file.
#'
#' @param obj h5Seurat filename or connection
#' @param ... additional parameters
#' @export
fetch_data <- function(obj, ...) {
  UseMethod("fetch_data")
}

#' @describeIn fetch_data Access data from filename.
#' @method fetch_data character
#' @export
fetch_data.character <- function(obj, ...) {
  hfile <- SeuratDisk::Connect("anterior1.h5Seurat")
  fetch_data(hfile, ...)
  hfile$close_all()
}

#' @describeIn fetch_data Access data from h5Seurat connection.
#' @method fetch_data h5Seurat
#' @export
fetch_data.h5Seurat <- function(obj, ...) {

}