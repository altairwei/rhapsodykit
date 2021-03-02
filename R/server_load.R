validate_manifest <- function(x) {
  stopifnot(is.list(x))

  if (is.null(x[[KEY_LIB_RESULTS]])) {
    stop(sprintf("'%s' in manifest.json is needed.", KEY_LIB_RESULTS))
  }
}

#' Load project manifest from root folder.
#'
#' @param base_folder Project results root folder.
#' @return A list converted from json file.
load_manifest <- function(base_folder) {
  if (!dir.exists(base_folder)) {
    stop("Directory provided does not exist: ", base_folder)
  }

  manifest_filename <- file.path(base_folder, "manifest.json")

  if (!file.exists(manifest_filename)) {
    stop("manifest.json was not found in ", base_folder)
  }

  # Read text first
  manifest <- rjson::fromJSON(
    readr::read_file(manifest_filename), simplify = FALSE)

  validate_manifest(manifest)

  manifest
}

#' Search Rhapsody WTA outputs for given folders
#'
#' @param raw_result_folders Result folders of Rhapsody WTA pipeline.
#' @param base_folder Prefix for result folders.
#' @return A list of results files for each sample.
search_data <- function(raw_result_folders, base_folder) {
  raw_result_folders <- file.path(base_folder, raw_result_folders)
  names(raw_result_folders) <- raw_result_folders
  libs <- raw_result_folders %>%
    lapply(function(lib_dir) {
      seurat_dir <- file.path(lib_dir, "SeuratAnalysis")
      list(
        Expression_Data = find_file(lib_dir, "*_Expression_Data.st"),
        Metrics_Summary = find_file(lib_dir, "*_Metrics_Summary.csv"),
        RSEC_MolsPerCell = find_file(lib_dir, "*_RSEC_MolsPerCell.csv"),
        RSEC_ReadsPerCell = find_file(lib_dir, "*_RSEC_ReadsPerCell.csv"),
        UMI_Adjusted_Stats = find_file(lib_dir, "*_UMI_Adjusted_Stats.csv"),
        Seurat_Object = find_file(seurat_dir, "Seurat_Object.rds")
      )
    })

  libs
}

load_expression_data <- function(filename) {
  df <- vroom::vroom(filename, delim = "\t", comment = "#",
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