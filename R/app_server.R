#' Main server logic
#'
#' @param input,output,session Arguments passed by shiny
main_server <- function(input, output, session) {
  data_folder <- shiny::getShinyOption("data_folder")

  manifest <- load_manifest(data_folder)

  rv <- shiny::reactiveValues()

  integrated_list <- search_integrated_data(
    manifest[[KEY_LIB_INTEGRATED]], data_folder)
  names(integrated_list) <- basename(names(integrated_list))

  library_list <- c(integrated_list)
  # A list of R object cache for each library
  cache <- shiny::reactiveValues(
    cell_embeddings = NULL,
    gene_expressions = NULL
  )

  #shiny::callModule(
  #    metrics_summary_page, "metrics_summary_page", rv, cache)
  shiny::callModule(
      expression_projection_page, "expression_projection_page",
      library_list, cache)
  shiny::callModule(
     subpopulation_heatmap_page, "subpopulation_heatmap_page",
     library_list, cache)
  # shiny::callModule(
  #     ensembl_homologs_page, "ensembl_homologs_page")

  waiter::waiter_hide()
}