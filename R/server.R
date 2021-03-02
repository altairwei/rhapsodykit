#' Main server logic
#'
#' @param input,output,session Arguments passed by shiny
main_server <- function(input, output, session) {
  data_folder <- shiny::getShinyOption("data_folder")

  manifest <- load_manifest(data_folder)

  rv <- shiny::reactiveValues()

  # A list of results files for each library.
  library_list <- search_data(manifest[[KEY_LIB_RESULTS]], data_folder)
  names(library_list) <- basename(names(library_list))

  # A list of R object cache for each library
  cache <- lapply(library_list, function(results_list) {
    purrr::imap(results_list, function(filename, key) {
      method_to_load <- switch(key,
        Expression_Data = load_expression_data,
        Seurat_Object = readRDS,
        function(x) NULL
      )

      shiny::reactive({
        waiter::waiter_show(html = tagList(
          waiter::spin_flower(),
          htmltools::h4(paste0("Reading ", basename(filename)))
        ))
        on.exit(waiter::waiter_hide(), add = TRUE)
        method_to_load(filename)
      })
    })
  })

  #shiny::callModule(
  #    metrics_summary_page, "metrics_summary_page", rv, cache)
  shiny::callModule(
      expression_projection_page, "expression_projection_page",
      library_list, cache)
}