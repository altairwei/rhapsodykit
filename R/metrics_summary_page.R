metrics_summary_page_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shiny::uiOutput(outputId = ns("cell_molecules_distribution"))
  )
}

metrics_summary_page <- function(input, output, session, rv) {
  # Render settings for scatter plot
  output$cell_molecules_distribution <- shiny::renderUI({
    shiny::req(rv$library_list)
    ns <- session$ns

    panel <- purrr::imap(rv$library_list, function(lib, name) {
      shiny::tabPanel(basename(name),
        shiny::plotOutput(outputId = ns("cell_mol_hist_plot")))
    })
    names(panel) <- NULL

    tabox_args <- list(width = 6)
    tabox_args <- c(tabox_args, panel)

    do.call(shinydashboard::tabBox, tabox_args)
  })
}