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
          choiceValues = c("umap", "tsne", "pca"),
          choiceNames = c("UMAP", "t-SNE", "PCA"),
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
    gene_list <- strsplit(stringr::str_trim(input$gene_list), "\n")[[1]]
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

    expression <- fetch_gene_expression(
      library_list[[library]]$Seurat_Disk,
      gene_list,
      cache
    )

    embedding <- fetch_cell_embeddings(
      library_list[[library]]$Seurat_Disk, cache)

    dim_names <- paste0(
      switch(input$reduction,
        tsne = "tSNE_", umap = "UMAP_", pca = "PC_"),
      1:2
    )

    data <- cbind(
      embedding[, c(dim_names, "ident")],
      expression
    )

    data <- tidyr::gather(
      data, tidyselect::any_of(paste0("rna_", gene_list)),
      key = "gene", value = "expr"
    )

    aggregated <- data %>%
      dplyr::group_by(gene, ident) %>%
      dplyr::summarise(value = fun(expr)) %>%
      # Convert data frame to matrix
      tidyr::pivot_wider(names_from = ident, values_from = value)

    mat <- data.matrix(aggregated[, -1])
    rownames(mat) <- stringr::str_remove(aggregated[["gene"]], "rna_")

    mat
  })

  ht_obj_cache <- shiny::reactiveVal(NULL)
  ht_pos_cache <- shiny::reactiveVal(NULL)

  sub_ht_obj_cache <- shiny::reactiveVal(NULL)
  sub_ht_pos_cache <- shiny::reactiveVal(NULL)

  # shiny::observeEvent(input$submit, {
  #   ht_obj_cache(NULL)
  #   ht_pos_cache(NULL)
  #   sub_ht_obj_cache(NULL)
  #   sub_ht_pos_cache(NULL)
  # })

  brushed_indexes <- shiny::reactive({
    ht <- ht_obj_cache()
    ht_pos <- ht_pos_cache()

    if (shiny::isTruthy(input$heatmap_brush) &&
          shiny::isTruthy(ht) && shiny::isTruthy(ht_pos)) {
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

  brushed_mat <- shiny::reactive({
    index_list <- brushed_indexes()
    mat <- dataset()
    if (shiny::isTruthy(index_list) && shiny::isTruthy(mat))
      mat[index_list$row_index, index_list$col_index, drop = FALSE]
    else
      NULL
  })

  clicked_gene <- shiny::reactiveVal(NULL)

  # Update clicked gene from main heatmap
  shiny::observe({
    ht <- ht_obj_cache()
    ht_pos <- ht_pos_cache()
    mat <- dataset()

    if (shiny::isTruthy(input$heatmap_click)
          && shiny::isTruthy(ht)
          && shiny::isTruthy(ht_pos)
          && shiny::isTruthy(mat)) {

      pos <- get_pos_from_click(input$heatmap_click)
      selection <- InteractiveComplexHeatmap::selectPosition(
        ht, pos,
        mark = FALSE,
        ht_pos = ht_pos,
        verbose = FALSE,
        calibrate = FALSE
      )

      row_index <- unique(unlist(selection$row_index))
      col_index <- unique(unlist(selection$column_index))

      sub <- mat[row_index, col_index, drop = FALSE]

      if (nrow(sub) != 0)
        clicked_gene(paste0("rna_", rownames(sub)))
      else
        clicked_gene(NULL)

      sub

    } else {
      NULL
    }

  })

  # Update clicked gene from main sub-heatmap
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

      clicked_gene(paste0("rna_", rownames(sub)))

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
    mtx <- dataset()

    shiny::req(mtx)

    ht <- plot_subpopulation_heatmap(
      mtx,
      nrow(mtx) < MAX_GENE_NUM,
      col = colorRampPalette(
        c("grey", "white", "yellow", "red", "dark red"))(256),
      heatmap_legend_param = list(
        title = "Expression Level",
        title_position = "leftcenter-rot",
        legend_height = grid::unit(6, "cm")
      )
    )

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
    subm <- brushed_mat()

    if (shiny::isTruthy(index_list)
          && shiny::isTruthy(ht_list)
          && shiny::isTruthy(subm)) {

      #message("Plot brushed heatmap...")
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
      #message("Finished plot brushed heatmap in ", elapsed, " secs")

      ht_selected
    } else {
      plot_placeholder("Select heatmap\n to see sub-heatmap.")
    }
  })

  clicked_data <- shiny::reactive({
    selected_gene <- clicked_gene()
    if (shiny::isTruthy(selected_gene)) {
      dim_names <- paste0(
        switch(input$reduction,
          tsne = "tSNE_", umap = "UMAP_", pca = "PC_"),
        1:2
      )
      cbind(
        shiny::isolate(cache$cell_embeddings)[, c(dim_names, "ident"), drop = FALSE],
        shiny::isolate(cache$gene_expressions)[, selected_gene, drop = FALSE]
      )
    } else {
      NULL
    }
  })

  output$projection <- shiny::renderPlot({
    data_to_plot <- clicked_data()
    gene <- clicked_gene()
    if (shiny::isTruthy(data_to_plot) && shiny::isTruthy(gene)) {
      reduction <- input$reduction
      feature_scatter(data_to_plot, gene, reduction = reduction)
    } else {
      plot_placeholder(
        "Click heatmap to \nselect a gene to projection.")
    }
  })

  output$vlnplot <- shiny::renderPlot({
    data_to_plot <- clicked_data()
    gene <- clicked_gene()
    if (shiny::isTruthy(data_to_plot) && shiny::isTruthy(gene)) {
      feature_violin(data_to_plot, gene)
    } else {
      plot_placeholder(
        "Click heatmap to \nselect a gene to projection.")
    }
  })

  output$selected_info <- DT::renderDataTable({
    mat <- brushed_mat()
    shiny::req(mat)
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