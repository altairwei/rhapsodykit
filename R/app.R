#' Run Rhapsody WTA Viewer application
#'
#' @param data_folder A folder containing data processing results。
#'
#' @export
run_app <- function(data_folder = ".", host = "127.0.0.1", port = NULL) {
  shiny::shinyOptions(
    data_folder = data_folder)
  opts <- list(
    host = host,
    port = port
  )
  app <- shiny::shinyApp(
    ui = create_main_ui(),
    server = main_server,
    options = opts
  )

  shiny::runApp(app)
}