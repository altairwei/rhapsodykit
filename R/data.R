pluck_endswith <- function(x, suffix) {
  purrr::detect(x, function(filename) {
    stringr::str_ends(filename, suffix)
  })
}

search_data <- function(root_folder) {
  fs::dir_ls(root_folder) %>% # Library
    purrr::map(function(lib_dir) {
      files <- fs::dir_ls(lib_dir, type = "file")
      list(
        Expression_Data = pluck_endswith(files, "_Expression_Data.st"),
        Metrics_Summary = pluck_endswith(files, "_Metrics_Summary.csv"),
        RSEC_MolsPerCell = pluck_endswith(files, "_RSEC_MolsPerCell.csv"),
        RSEC_ReadsPerCell = pluck_endswith(files, "_RSEC_ReadsPerCell.csv"),
        UMI_Adjusted_Stats = pluck_endswith(files, "_UMI_Adjusted_Stats.csv")
      )
    })
}

load_expression_data <- function(filename) {
  df <- readr::read_tsv(filename, comment = "#",
    col_types = readr::cols(
      Cell_Index = readr::col_factor(),
      Gene = readr::col_character(),
      RSEC_Reads = readr::col_double(),
      Raw_Molecules = readr::col_double(),
      RSEC_Adjusted_Molecules = readr::col_double()
    )
  )

  df
}

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

read_rhapsody_wta <- function(base_dir, use_mtx = FALSE) {
  expr_matrix <- NULL

  if (isTRUE(use_mtx)) {
    expr_matrix <- read_mtx(base_dir)
  } else {
    expr_matrix <- read_raw_csv(base_dir)
  }

  expr_matrix
}