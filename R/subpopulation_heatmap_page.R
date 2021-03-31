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
      shinydashboard::box(
        width = NULL,
        shinycssloaders::withSpinner(
          shiny::plotOutput(outputId = ns("projection"))
        )
      )
    ),
    shiny::column(width = 6,
      shiny::uiOutput(outputId = ns("heatmap_container"))
    ),
    shiny::column(width = 3,
      shinydashboard::box(
        width = NULL,
        title = "Selected Data",
        DT::dataTableOutput(
          outputId = ns("selected_info")
        )
      )
    ),
  )
}

subpopulation_heatmap_page <-  function(
  input, output, session, library_list, cache) {

  shiny::updateSelectInput(
    session,
    "select_sample",
    choices = basename(names(library_list))
  )

  ####################
  # User inputs
  ####################

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

  ####################
  # Cache Data
  ####################

  dataset <- shiny::reactive({
    gene_list <- gene_queries()
    fun <- shiny::isolate(stat_fun())

    library <- shiny::isolate(input$select_sample)
    obj <- cache[[library]]$Seurat_Object()

    aggregated <- aggregate_by_ident(obj, gene_list, fun)

    # Convert data frame to matrix
    aggregated <- aggregated %>%
      tidyr::pivot_wider(names_from = ident, values_from = value)
    mat <- data.matrix(aggregated[, -1])
    rownames(mat) <- aggregated[["gene"]]

    mat
  })

  ht_obj_cache <- shiny::reactiveVal(NULL)
  ht_pos_cache <- shiny::reactiveVal(NULL)

  clicked_mat <- shiny::reactive({

    ht <- ht_obj_cache()
    ht_pos <- ht_pos_cache()

    if (shiny::isTruthy(input$heatmap_click) &&
      shiny::isTruthy(ht) && shiny::isTruthy(ht_pos)) {
      pos <- get_pos_from_click(input$heatmap_click)
      selection <- InteractiveComplexHeatmap::selectPosition(
        ht, pos,
        mark = FALSE,
        ht_pos = ht_pos_cache(),
        verbose = FALSE,
        calibrate = FALSE
      )

      row_index <- unique(unlist(selection$row_index))
      col_index <- unique(unlist(selection$column_index))

      mat <- dataset()

      mat[row_index, col_index, drop = FALSE]

    } else {
      NULL
    }

  })

  brushed_mat <- shiny::eventReactive(input$heatmap_brush, {
    shiny::req(ht_obj_cache(), ht_pos_cache())
    ht <- ht_obj_cache()

    lt <- get_pos_from_brush(input$heatmap_brush)
    pos1 <- lt[[1]]
    pos2 <- lt[[2]]

    selection <- InteractiveComplexHeatmap::selectArea(
      ht,
      mark = FALSE,
      pos1 = pos1,
      pos2 = pos2,
      verbose = FALSE,
      ht_pos = ht_pos_cache(),
      include_annotation = TRUE,
      calibrate = FALSE
    )

    row_index <- unique(unlist(selection$row_index))
    col_index <- unique(unlist(selection$column_index))

    mat <- dataset()

    mat[row_index, col_index, drop = FALSE]
  })

  ####################
  # Response
  ####################

  #TODO: 点击和框选事件在更新热图后会失效
  #TODO: 允许拷贝整个表格，而非特定页

  output$heatmap_container <- shiny::renderUI({
    ns <- session$ns
    gene_list <- gene_queries()

    height_per_gene <- 14
    height <- ifelse(
      length(gene_list) > MAX_GENE_NUM,
      50 + height_per_gene * MAX_GENE_NUM,
      50 + height_per_gene * length(gene_list))

    shinydashboard::box(
      width = NULL, title = "Heatmap",
      shinycssloaders::withSpinner(
        shiny::plotOutput(
          outputId = ns("heatmap"),
          height = sprintf("%ipx", height),
          click = ns("heatmap_click"),
          brush = ns("heatmap_brush"))
      )
    )

  })

  output$heatmap <- shiny::renderPlot({
    gene_list <- gene_queries()

    ht <- plot_subpopulation_heatmap(
      dataset(), length(gene_list) < MAX_GENE_NUM)
    ht <- ComplexHeatmap::draw(ht)
    ht_obj_cache(ht)

    ht_pos <- InteractiveComplexHeatmap::htPositionsOnDevice(
      ht, include_annotation = TRUE, calibrate = FALSE)
    ht_pos_cache(ht_pos)

    ht
  })

  output$projection <- shiny::renderPlot({
    mat <- clicked_mat()

    if (shiny::isTruthy(mat)) {
      selected_gene <- rownames(mat)

      library <- shiny::isolate(input$select_sample)
      reduction <- input$reduction
      obj <- cache[[library]]$Seurat_Object()
      Seurat::DefaultAssay(obj) <- "RNA"
      Seurat::FeaturePlot(obj, features = selected_gene, reduction = reduction)
    } else {
      grid::grid.text(
        "Click heatmap to \nselect a gene to projection.",
        0.5, 0.5, gp = grid::gpar(fontsize = 20))
    }
  })

  output$selected_info <- DT::renderDataTable({
    mat <- brushed_mat()
    mat <- round(mat, 2)
    DT::datatable(
      mat,
      rownames = TRUE,
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
              columns = 0,
              modifier = list(
                page = "all"
              )
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
  })
}