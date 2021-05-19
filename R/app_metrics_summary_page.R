metrics_summary_page_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shiny::uiOutput(outputId = ns("cell_molecules_distribution"))
  )
}

metrics_summary_page <- function(input, output, session, rv, cache) {

  output$cell_molecules_distribution <- shiny::renderUI({
    shiny::req(rv$library_list)
    ns <- session$ns

    panel <- purrr::imap(rv$library_list, function(lib, name) {
      shiny::tabPanel(basename(name),
        shiny::plotOutput(
          outputId = ns(
            paste("cell_mol_hist_plot", basename(name), sep = "_"))))
    })
    names(panel) <- NULL

    tabox_args <- list(width = 6)
    tabox_args <- c(tabox_args, panel)

    do.call(shinydashboard::tabBox, tabox_args)
  })

  shiny::observe({
    shiny::req(rv$library_list)
    ns <- session$ns

    for (libpath in names(rv$library_list)) {
      # We `local()` to capture `libpath` for `shiny::renderPlot` expressions,
      #   otherwise `libpath` will be a fixed value when loop ended.
      local({
        plot_ui_id <- paste("cell_mol_hist_plot", basename(libpath), sep = "_")
        lib_id <- libpath
        output[[plot_ui_id]] <- shiny::renderPlot({
          cache$Expression_Data()[[lib_id]] %>%
            dplyr::group_by(Cell_Index) %>%
            dplyr::summarise(Molecules = sum(RSEC_Adjusted_Molecules)) %>%
            ggplot2::ggplot(ggplot2::aes(Molecules)) +
            ggplot2::geom_histogram(bins = 100)
        })
      })
    }
  })
}