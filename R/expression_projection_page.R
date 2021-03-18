expression_projection_page_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shiny::column(width = 3,
      shinydashboard::box(width = NULL,
        shiny::selectInput(
          inputId = ns("select_sample"),
          label = "Select Sample:",
          choices = c("Please select a sample..." = "")
        ),
        shiny::radioButtons(
          inputId = ns("reduction"),
          label = "Reduction",
          choices = c("tsne", "umap", "pca"),
          inline = TRUE
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
      ),
      shinycssloaders::withSpinner(
        shiny::plotOutput(
          outputId = ns("dimplot")
        )
      )
    ),
    shiny::column(width = 9,
      shiny::uiOutput(outputId = ns("display"))
    )
  )
}

expression_projection_page <- function(
  input, output, session, library_list, cache) {
  shiny::updateSelectInput(
    session,
    "select_sample",
    choices = basename(names(library_list)),
  )

  gene_queries <- shiny::eventReactive(input$submit, {
    gene_list <- strsplit(input$gene_list, "\n")[[1]]
    gene_list
  })

  output$dimplot <- shiny::renderPlot({
    shiny::req(input$submit)
    shiny::req(input$select_sample)

    library <- input$select_sample
    reduction <- input$reduction
    obj <- cache[[library]]$Seurat_Object()

    Seurat::DefaultAssay(obj) <- "RNA"
    Seurat::DimPlot(obj, reduction = reduction)
  })

  # Dynamically render plotOutput that user need.
  output$display <- shiny::renderUI({
    ns <- session$ns
    gene_list <- gene_queries()

    lapply(gene_list, function(gene) {
      shinydashboard::box(width = NULL, title = gene,
        shiny::column(6,
          shinycssloaders::withSpinner(
            shiny::plotOutput(outputId = ns(paste0("scatter-", gene)))
          )
        ),
        shiny::column(6,
          shinycssloaders::withSpinner(
            shiny::plotOutput(outputId = ns(paste0("violin-", gene)))
          )
        )
      )
    })
  })

  last_queries <- NULL
  shiny::observe({
    # Remove last expired output observer.
    if (!is.null(last_queries)) {
      for (old_id in last_queries) {
        output[[paste0("scatter-", old_id)]] <- NULL
        output[[paste0("violin-", old_id)]] <- NULL
      }
    }

    plot_output_id_list <- gene_queries()
    last_queries <<- plot_output_id_list

    library <- shiny::isolate(input$select_sample)
    reduction <- shiny::isolate(input$reduction)
    obj <- cache[[library]]$Seurat_Object()

    invisible(lapply(plot_output_id_list, function(gene_id) {
      output[[paste0("scatter-", gene_id)]] <- shiny::renderPlot({
        Seurat::FeaturePlot(obj, features = gene_id, reduction = reduction)
      })
      output[[paste0("violin-", gene_id)]] <- shiny::renderPlot({
        Seurat::VlnPlot(obj, features = gene_id)
      })
    }))
  })
}