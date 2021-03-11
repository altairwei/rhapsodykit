ensembl_homologs_page_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shiny::column(width = 3,
      shinydashboard::box(width = NULL,
        id = ns("controller"),
        shiny::selectInput(
          inputId = ns("host"),
          label = "Select host:",
          choices = c(
            "CHOOSE HOST" = "",
            "Ensembl Plants" = "plants.ensembl.org"
            )
        ),
        shiny::selectInput(
          inputId = ns("mart"),
          label = "Select mart:",
          choices = c("CHOOSE MART" = "")
        ),
        shiny::selectInput(
          inputId = ns("dataset"),
          label = "Select dataset:",
          choices = c("CHOOSE DATASET" = "")
        ),
        shiny::selectInput(
          inputId = ns("target"),
          label = "Select target species:",
          choices = c("CHOOSE TARGET" = "")
        ),
        shiny::textInput(
          inputId = ns("query_gene"),
          label = "Query gene:"
        ),
        shiny::actionButton(
          inputId = ns("submit"),
          label = "Submit"
        )
      )
    ),
    shiny::column(width = 9,
      shinycssloaders::withSpinner(
        DT::dataTableOutput(
          outputId = ns("homolog_table")
        )
      )
    )
  )
}

ensembl_homologs_page <- function(input, output, session) {
    ns <- session$ns

    wtr <- waiter::Waiter$new(
      id = ns("controller"),
      html = waiter::spin_throbber(),
      color = "rgba(255,255,255,.5)")

    # Update available marts
    shiny::observeEvent(input$host, {
      shiny::req(input$host)

      wtr$show()
      on.exit(wtr$hide(), add = TRUE)

      all_marts <- biomaRt::listMarts(host = input$host)

      choices <- all_marts$biomart
      names(choices) <- all_marts$version

      shiny::updateSelectInput(
        session, "mart", choices = choices)
    })

    # Update available datasets
    shiny::observeEvent(input$mart, {
      shiny::req(input$mart)

      wtr$show()
      on.exit(wtr$hide(), add = TRUE)

      ensembl <- biomaRt::useMart(input$mart, host = input$host)
      all_datasets <- biomaRt::listDatasets(ensembl)

      choices <- all_datasets$dataset
      names(choices) <- all_datasets$description

      shiny::updateSelectInput(
        session, "dataset", choices = choices)
    })

    # Update available target species
    shiny::observeEvent(input$dataset, {
      shiny::req(input$dataset)

      wtr$show()
      on.exit(wtr$hide(), add = TRUE)

      dataset_conn <- biomaRt::useMart(
        input$mart, dataset = input$dataset, host = input$host)

      species <- biomaRt::listAttributes(dataset_conn) %>%
        dplyr::filter(
          page == "homologs",
          stringr::str_ends(name, "_homolog_ensembl_gene")) %>%
        dplyr::transmute(
          name = stringr::str_remove(name, "_homolog_ensembl_gene"),
          description = stringr::str_remove(description, " gene stable ID"))

      choices <- species$name
      names(choices) <- species$description

      shiny::updateSelectInput(
        session, "target", choices = choices)
    })

    output$homolog_table <- DT::renderDataTable(
      {
        shiny::req(input$submit)

        dataset_conn <- biomaRt::useMart(
          shiny::isolate(input$mart),
          dataset = shiny::isolate(input$dataset),
          host = shiny::isolate(input$host)
        )

        attrs_to_retrive <- c(
          "_homolog_ensembl_gene",
          "_homolog_orthology_type",
          "_homolog_perc_id",
          "_homolog_perc_id_r1",
          "_homolog_wga_coverage",
          "_homolog_orthology_confidence"
        )

        attrs_to_retrive <- paste0(
          shiny::isolate(input$target), attrs_to_retrive)

        names(attrs_to_retrive) <- c(
            "Gene ID", "Type", "Target %id", "Query %id",
            "WGA Coverage", "Confidence"
        )

        all_homologs <- biomaRt::getBM(
          attributes = attrs_to_retrive,
          filters = "ensembl_gene_id",
          values = shiny::isolate(input$query_gene),
          mart = dataset_conn
        )

        all_homologs <- all_homologs %>%
          dplyr::rename_with(~ names(attrs_to_retrive))

        DT::datatable(
          all_homologs,
          selection = list(target = "column"),
          extensions = "Buttons",
          options = list(
            pageLength = 50,
            scrollX = TRUE,
            dom = "Bfrtip",
            buttons = list(
              "pageLength",
              list(
                extend = "copyHtml5",
                text = "Copy Gene ID",
                title = NULL,
                header = FALSE,
                exportOptions = list(
                    columns = 1
                )
              ),
              list(
                extend = "copyHtml5",
                text = "Copy Table",
                title = NULL
              ),
              "csv",
              "print"
            )
          )
        )
      },
      server = TRUE
    )
}