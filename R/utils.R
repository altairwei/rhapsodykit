pluck_endswith <- function(x, suffix) {
  purrr::detect(x, function(filename) {
    stringr::str_ends(filename, suffix)
  })
}

find_file <- function(base_folder, glob) {
  if (!dir.exists(base_folder)) {
    stop("Directory provided does not exist: ", base_folder)
  }

  filename <- Sys.glob(file.path(base_folder, glob))

  if (length(filename) < 1) {
    stop(sprintf("can not find `%s` file", glob))
  } else if (length(filename) > 1) {
    stop(sprintf("more than one `%s` file was found.", glob))
  }

  filename
}