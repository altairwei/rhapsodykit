#' Main server logic
#'
#' @param input,output,session Arguments passed by shiny
main_server <- function(input, output, session) {
  data_folder <- shiny::getShinyOption("data_folder")
  rv <- shiny::reactiveValues(
    library_list = search_data(data_folder),
    cache = list()
  )

  shiny::callModule(
      metrics_summary_page, "metrics_summary_page", rv)

}