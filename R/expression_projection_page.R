expression_projection_page_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shiny::column(width = 3,
      shinydashboard::box(width = NULL,
        shiny::selectInput(
          inputId = ns("select_sample"),
          label = "Select Sample:",
          choices = c("Please select a sample..." = "None")
        ),
        shiny::textAreaInput(
          inputId = ns("gene_list"),
          label = "Genes to Query:",
          height = "200px",
          placeholder = paste0("Please enter the list of genes",
            " to be queried, one gene per line.")
        ),
        shiny::actionButton(
          inputId = ns("submit"),
          label = "Submit"
        )
      )
    ),
    shiny::column(width = 9,
      shiny::uiOutput(outputId = ns("display"))
    )
  )
}

expression_projection_page <- function(input, output, session, rv, cache) {
  shiny::observe({
    shiny::req(rv$library_list)
    shiny::updateSelectInput(
      session,
      "select_sample",
      choices = basename(names(rv$library_list)),
    )
  })

  gene_queries <- shiny::eventReactive(input$submit, {
    gene_list <- strsplit(input$gene_list, "\n")[[1]]
    gene_list
  })

  output$display <- shiny::renderUI({
    ns <- session$ns
    gene_list <- gene_queries()

    lapply(gene_list, function(gene) {
      shinydashboard::box(width = NULL, title = gene,
        shiny::column(6,
          shiny::plotOutput(outputId = ns(paste0("scatter-", gene)))
        ),
        shiny::column(6,
          shiny::plotOutput(outputId = ns(paste0("violin-", gene)))
        )
      )
    })
  })

  last_queries <- NULL
  shiny::observe({
    if (!is.null(last_queries)) {
      for (old_id in last_queries) {
        output[[paste0("scatter-", old_id)]] <- NULL
        output[[paste0("violin-", old_id)]] <- NULL
      }
    }

    plot_output_id_list <- gene_queries()
    last_queries <<- plot_output_id_list

    for (id in plot_output_id_list) {
      data <- data.frame(
        x = rnorm(nchar(id), 1, 2),
        y = rnorm(nchar(id), 1, 2)
      )
      output[[paste0("scatter-", id)]] <- shiny::renderPlot({
        ggplot2::ggplot(data, ggplot2::aes(x, y)) +
          ggplot2::geom_point()
      })
      output[[paste0("violin-", id)]] <- shiny::renderPlot({
        ggplot2::ggplot(data, ggplot2::aes(y, x)) +
          ggplot2::geom_point()
      })
    }
  })
}