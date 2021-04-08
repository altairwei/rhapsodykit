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
        ),
        shinycssloaders::withSpinner(
          shiny::plotOutput(outputId = ns("vlnplot"))
        )
      )
    ),
    shiny::column(width = 6,
      shiny::uiOutput(outputId = ns("heatmap_container"))
    ),
    shiny::column(width = 3,
      shinydashboard::box(
        width = NULL,
        title = "Sub-heatmap",
        shinycssloaders::withSpinner(
          shiny::plotOutput(
            outputId = ns("sub_heatmap"),
            click = ns("sub_heatmap_click"))
        )
      ),
      shinydashboard::box(
        width = NULL,
        title = "Selected Data",
        shinycssloaders::withSpinner(
          DT::dataTableOutput(
            outputId = ns("selected_info")
          )
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

  sub_ht_obj_cache <- shiny::reactiveVal(NULL)
  sub_ht_pos_cache <- shiny::reactiveVal(NULL)

  brushed_indexes <- shiny::reactive({
    ht <- ht_obj_cache()
    ht_pos <- ht_pos_cache()

    if (shiny::isTruthy(input$heatmap_brush) &&
          !is.null(ht) && !is.null(ht_pos)) {
      lt <- get_pos_from_brush(input$heatmap_brush)
      pos1 <- lt[[1]]
      pos2 <- lt[[2]]

      selection <- InteractiveComplexHeatmap::selectArea(
        ht,
        mark = FALSE,
        pos1 = pos1,
        pos2 = pos2,
        verbose = FALSE,
        ht_pos = ht_pos,
        include_annotation = TRUE,
        calibrate = FALSE
      )

      row_index <- unique(unlist(selection$row_index))
      col_index <- unique(unlist(selection$column_index))

      list(
        row_index = row_index,
        col_index = col_index,
        heatmap = selection$heatmap
      )

    } else {
      NULL
    }
  })

  brushed_mat <- shiny::eventReactive(input$heatmap_brush, {
    index_list <- brushed_indexes()

    shiny::req(index_list)

    mat <- dataset()

    mat[index_list$row_index, index_list$col_index, drop = FALSE]
  })

  clicked_gene <- shiny::reactiveVal(NULL)

  shiny::observe({
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

      sub <- mat[row_index, col_index, drop = FALSE]

      clicked_gene(rownames(sub))

      sub

    } else {
      NULL
    }

  })

  shiny::observe({
    ht <- sub_ht_obj_cache()
    ht_pos <- sub_ht_pos_cache()

    if (shiny::isTruthy(input$sub_heatmap_click) &&
      shiny::isTruthy(ht) && shiny::isTruthy(ht_pos)) {
      pos <- get_pos_from_click(input$sub_heatmap_click)
      selection <- InteractiveComplexHeatmap::selectPosition(
        ht, pos,
        mark = FALSE,
        ht_pos = ht_pos,
        verbose = FALSE,
        calibrate = FALSE
      )

      row_index <- unique(unlist(selection$row_index))
      col_index <- unique(unlist(selection$column_index))

      mat <- brushed_mat()

      sub <- mat[row_index, col_index, drop = FALSE]

      clicked_gene(rownames(sub))

      sub
    } else {
      NULL
    }
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

  ht_raw_obj <- shiny::reactive({
    gene_list <- gene_queries()

    ht <- plot_subpopulation_heatmap(
      dataset(), length(gene_list) < MAX_GENE_NUM)

    ht
  })

  output$heatmap <- shiny::renderPlot({
    ht <- ht_raw_obj()
    ht <- ComplexHeatmap::draw(ht)
    ht_obj_cache(ht)

    ht_pos <- InteractiveComplexHeatmap::htPositionsOnDevice(
      ht, include_annotation = TRUE, calibrate = FALSE)
    ht_pos_cache(ht_pos)

    ht
  })

  output$sub_heatmap <- shiny::renderPlot({
    index_list <- brushed_indexes()
    ht_list <- ht_obj_cache()

    if (shiny::isTruthy(index_list) && !is.null(ht_list)) {
      subm <- brushed_mat()

      message("Plot brushed heatmap...")
      start_time <- Sys.time()

      ri <- index_list$row_index
      ci <- index_list$col_index

      ht_current_full <- ht_list@ht_list[[index_list$heatmap[[1]]]]

      row_labels <- ht_current_full@row_names_param$labels
      if (!is.null(row_labels)) {
        row_labels <- row_labels[ri]
      }
      column_labels <- ht_current_full@column_names_param$labels
      if (!is.null(column_labels)) {
        column_labels <- column_labels[ci]
      }

      ht_selected <- ComplexHeatmap::Heatmap(
        subm,
        rect_gp = ht_current_full@matrix_param$gp,
        col = ht_current_full@matrix_color_mapping,
        show_heatmap_legend = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_title = NULL,
        column_title = NULL,
        border = ht_current_full@matrix_param$border,
        row_labels = row_labels,
        column_labels = column_labels,
        show_row_names = length(ri) < 50,
        row_names_side = ht_current_full@row_names_param$side,
        show_column_names = TRUE,
        column_names_side = ht_current_full@column_names_param$side
      )

      ht_selected <- ComplexHeatmap::draw(ht_selected)
      sub_ht_obj_cache(ht_selected)

      ht_selected_pos <- InteractiveComplexHeatmap::htPositionsOnDevice(
        ht_selected, include_annotation = TRUE, calibrate = FALSE)
      sub_ht_pos_cache(ht_selected_pos)

      end_time <- Sys.time()
      elapsed <- end_time - start_time
      message("Finished plot brushed heatmap in ", elapsed, " secs")

      ht_selected
    } else {
      plot_placeholder("Select heatmap\n to see sub-heatmap.")
    }
  })

  output$projection <- shiny::renderPlot({
    selected_gene <- clicked_gene()
    if (!is.null(selected_gene)) {
      library <- shiny::isolate(input$select_sample)
      reduction <- input$reduction
      obj <- cache[[library]]$Seurat_Object()
      Seurat::DefaultAssay(obj) <- "RNA"
      Seurat::FeaturePlot(obj, features = selected_gene, reduction = reduction)
    } else {
      plot_placeholder(
        "Click heatmap to \nselect a gene to projection.")
    }
  })

  output$vlnplot <- shiny::renderPlot({
    selected_gene <- clicked_gene()
    if (!is.null(selected_gene)) {
      library <- shiny::isolate(input$select_sample)
      reduction <- input$reduction
      obj <- cache[[library]]$Seurat_Object()
      Seurat::DefaultAssay(obj) <- "RNA"
      Seurat::VlnPlot(obj, features = selected_gene)
    } else {
      plot_placeholder(
        "Click heatmap to \nselect a gene to projection.")
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