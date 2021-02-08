validate_manifest <- function(x) {
  stopifnot(is.list(x))

  if (is.null(x[[KEY_LIB_RESULTS]])) {
    stop(sprintf("'%s' in manifest.json is needed.", KEY_LIB_RESULTS))
  }
}

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
  libs <- raw_result_folders %>% # Library
    lapply(function(lib_dir) {
      files <- fs::dir_ls(lib_dir, type = "file")
      list(
        Expression_Data = pluck_endswith(files, "_Expression_Data.st"),
        Metrics_Summary = pluck_endswith(files, "_Metrics_Summary.csv"),
        RSEC_MolsPerCell = pluck_endswith(files, "_RSEC_MolsPerCell.csv"),
        RSEC_ReadsPerCell = pluck_endswith(files, "_RSEC_ReadsPerCell.csv"),
        UMI_Adjusted_Stats = pluck_endswith(files, "_UMI_Adjusted_Stats.csv")
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