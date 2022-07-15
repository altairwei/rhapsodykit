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
            "Ensembl Plants" = "https://plants.ensembl.org"
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
        shiny::textAreaInput(
          inputId = ns("query_gene"),
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
      choices <- c("CHOOSE MART" = "", choices)

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
      choices <- c("CHOOSE DATASET" = "", choices)

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
      choices <- c("CHOOSE TARGET" = "", choices)

      shiny::updateSelectInput(
        session, "target", choices = choices)
    })

    gene_queries <- shiny::eventReactive(input$submit, {
      gene_list <- strsplit(stringr::str_trim(input$query_gene), "\n")[[1]]
      unique(gene_list)
    })

    output$homolog_table <- DT::renderDataTable(
      {
        shiny::req(input$submit)

        dataset_conn <- biomaRt::useMart(
          shiny::isolate(input$mart),
          dataset = shiny::isolate(input$dataset),
          host = shiny::isolate(input$host)
        )

        available_attrs <- biomaRt::listAttributes(dataset_conn)$name

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

        attrs_to_retrive <- c("ensembl_gene_id", attrs_to_retrive)
        names(attrs_to_retrive) <- c(
            "Query ID", "Gene ID", "Type",
            "Target %id", "Query %id",
            "WGA Coverage", "Confidence"
        )

        all_homologs <- biomaRt::getBM(
          attributes = intersect(attrs_to_retrive, available_attrs),
          filters = "ensembl_gene_id",
          values = gene_queries(),
          mart = dataset_conn
        )

        if (nrow(all_homologs) > 0) {
          all_homologs[setdiff(attrs_to_retrive, available_attrs)] <- "NA"
        } else {
          all_homologs[setdiff(attrs_to_retrive, available_attrs)] <- logical(0)
        }

        all_homologs <- all_homologs[attrs_to_retrive]

        all_homologs <- all_homologs %>%
          dplyr::rename_with(~ names(attrs_to_retrive)) %>%
          dplyr::filter(`Gene ID` != "")

        DT::datatable(
          all_homologs,
          rownames = FALSE,
          extensions = "Buttons",
          options = list(
            pageLength = 50,
            scrollX = TRUE,
            dom = "Bfrtip",
            buttons = list(
              "pageLength",
              list(
                extend = "copyHtml5",
                text = "Copy Target Gene ID",
                title = NULL,
                header = FALSE,
                exportOptions = list(
                  columns = 1,
                  modifier = list(page = "all")
                )
              ),
              list(
                extend = "copyHtml5",
                text = "Copy Table",
                title = NULL,
                exportOptions = list(
                  modifier = list(page = "all")
                )
              ),
              list(
                extend = "csv",
                title = NULL,
                exportOptions = list(
                  modifier = list(page = "all")
                )
              )
            )
          )
        )
      },
      # If server=FALSE then the buttons will export all data in the table,
      # while if server=TRUE they will only export visible data.
      server = FALSE
    )
}