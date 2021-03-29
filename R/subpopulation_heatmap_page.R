MAX_GENE_NUM <- 150

subpopulation_heatmap_page_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shiny::column(width = 3,
      shinydashboard::box(width = NULL, title = "Subpopulation Heatmap",
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
        shiny::radioButtons(
          inputId = ns("value_type"),
          label = "Value Type",
          choices = c("Mean", "Median", "SD"),
          inline = TRUE
        ),
        shiny::textAreaInput(
          inputId = ns("gene_list"),
          label = "Genes to Query:",
          height = "200px",
          placeholder = sprintf(
            paste0("Please enter the list of genes",
            " to be queried, one gene per line. If the number ",
            "of genes is greater than %i, then the gene name ",
            "will not be displayed."), MAX_GENE_NUM)
        ),
        shiny::actionButton(
          inputId = ns("submit"),
          label = "Submit"
        )
      ),
      shiny::plotOutput(outputId = ns("projection"))
    ),
    shiny::column(width = 9,
      shiny::uiOutput(outputId = ns("heatmap_container"))
    )
  )
}

subpopulation_heatmap_page <-  function(
  input, output, session, library_list, cache) {
  shiny::updateSelectInput(
    session,
    "select_sample",
    choices = basename(names(library_list))
  )

  gene_queries <- shiny::eventReactive(input$submit, {
    gene_list <- strsplit(input$gene_list, "\n")[[1]]
    unique(gene_list)
  })

  stat_fun <- shiny::reactive({
    switch(
      input$value_type,
      Mean = mean,
      Median = median,
      SD = sd,
      mean
    )
  })

  dataset <- shiny::reactive({
    gene_list <- gene_queries()
    fun <- shiny::isolate(stat_fun())

    library <- shiny::isolate(input$select_sample)
    obj <- cache[[library]]$Seurat_Object()

    aggregate_by_ident(obj, gene_list, fun)
  })

  output$heatmap_container <- shiny::renderUI({
    ns <- session$ns
    gene_list <- gene_queries()

    height_per_gene <- 14
    height <- ifelse(
      length(gene_list) > MAX_GENE_NUM,
      height_per_gene * MAX_GENE_NUM,
      height_per_gene * length(gene_list))

    shinycssloaders::withSpinner(
      shiny::plotOutput(
        outputId = ns("heatmap"),
        height = sprintf("%ipx", height),
        click = ns("heatmap_click"))
    )
  })

  output$heatmap <- shiny::renderPlot({
    shiny::req(input$submit)

    gene_list <- gene_queries()

    plot_subpopulation_heatmap(
      dataset(), length(gene_list) < MAX_GENE_NUM)
  })

  output$projection <- shiny::renderPlot({
    shiny::req(input$heatmap_click)

    df <- dataset()
    clicked_df <- shiny::nearPoints(
      df, input$heatmap_click, allRows = TRUE, addDist = TRUE)

    selected_gene <- clicked_df %>%
      dplyr::arrange(dist_) %>%
      .[[1, "gene"]]

    library <- shiny::isolate(input$select_sample)
    reduction <- input$reduction
    obj <- cache[[library]]$Seurat_Object()
    Seurat::DefaultAssay(obj) <- "RNA"
    Seurat::FeaturePlot(obj, features = selected_gene, reduction = reduction)
  })
}