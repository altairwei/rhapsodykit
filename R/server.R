#' Main server logic
#'
#' @param input,output,session Arguments passed by shiny
main_server <- function(input, output, session) {
  data_folder <- shiny::getShinyOption("data_folder")
  rv <- shiny::reactiveValues(
    library_list = search_data(data_folder)
  )

  cache <- list(
    Expression_Data = shiny::reactive({
      data <- NULL
      shiny::withProgress(message = "Loading Expression Data", value = 0, {
        n <- length(rv$library_list)
        data <- purrr::imap(rv$library_list, function(lib, name) {
          shiny::incProgress(1 / n, detail = basename(name))
          load_expression_data(lib[["Expression_Data"]])
        })
      })
      data
    })
  )

  shiny::callModule(
      metrics_summary_page, "metrics_summary_page", rv, cache)

}