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