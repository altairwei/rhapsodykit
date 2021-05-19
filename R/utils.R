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

#' Mapping key according to a given dict.
#'
#' @param x A character vector to convert.
#' @param dict A named character defined mapping rule.
#'
#' @return A character vector.
key_mapping <- function(x, dict) {
  stopifnot(is.character(dict), !is.null(names(dict)))

  available_keys <- names(dict)
  ret <- lapply(x, function(key) {
    if (key %in% available_keys) {
      dict[[key]]
    } else {
      key
    }
  })

  unlist(ret)
}
