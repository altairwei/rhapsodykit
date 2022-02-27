expression_projection_page_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shiny::column(width = 3,
      shinydashboard::box(width = NULL, title = "Expression Projection",
        shiny::selectInput(
          inputId = ns("select_sample"),
          label = "Select Sample:",
          choices = c("Please select a sample..." = "")
        ),
        shiny::radioButtons(
          inputId = ns("reduction"),
          label = "Dimension Reduction",
          choiceValues = c("umap", "tsne", "pca"),
          choiceNames = c("UMAP", "t-SNE", "PCA"),
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

fetch_gene_expression <- function(h5seurat_file, genes, cache) {
  if (!is.null(shiny::isolate(cache$gene_expressions))) {
    genes_to_query <- setdiff(
      paste0("rna_", genes),
      names(shiny::isolate(cache$gene_expressions))
    )
  } else {
    genes_to_query <- paste0("rna_", genes)
  }

  if (length(genes_to_query) > 0) {
    waiter::waiter_show(html = tagList(
      waiter::spin_flower(),
      htmltools::h4("Loading Gene Expressions..."),
      htmltools::p(paste(genes_to_query, collapse = " "))
    ))
    on.exit(waiter::waiter_hide(), add = TRUE)

    hfile <- SeuratDisk::Connect(h5seurat_file)
    on.exit(hfile$close_all(), add = TRUE)
    expressions <- SeuratDisk::FetchCellData(hfile, vars = genes_to_query)

    if (!is.null(shiny::isolate(cache$gene_expressions)))
      cache$gene_expressions[genes_to_query] <- expressions
    else
      cache$gene_expressions <- expressions
  }

  exprs <- shiny::isolate(cache$gene_expressions)
  genes_found <- intersect(paste0("rna_", genes), names(exprs))
  exprs[, genes_found]
}

fetch_cell_embeddings <- function(h5seurat_file, cache) {
  if (is.null(shiny::isolate(cache$cell_embeddings))) {
    waiter::waiter_show(html = tagList(
      waiter::spin_flower(),
      htmltools::h4("Loading Cell Embeddings...")
    ))
    on.exit(waiter::waiter_hide(), add = TRUE)

    hfile <- SeuratDisk::Connect(h5seurat_file)
    on.exit(hfile$close_all(), add = TRUE)

    embeddings <- SeuratDisk::FetchCellData(
      object = hfile,
      vars = c(
        "PC_1", "PC_2",
        "tSNE_1", "tSNE_2",
        "UMAP_1", "UMAP_2",
        "ident"
      )
    )

    cache$cell_embeddings <- embeddings
    embeddings
  } else {
    shiny::isolate(cache$cell_embeddings)
  }
}

cluster_plot <- function(data, reduction) {
  dim_names <- paste0(
    switch(reduction,
      tsne = "tSNE_", umap = "UMAP_", pca = "PC_"),
    1:2
  )

  ggplot2::ggplot(
      data = data,
      mapping = ggplot2::aes_string(
        x = dim_names[1], y = dim_names[2], color = "ident")) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_color_manual(
      values = Seurat::DiscretePalette(length(levels(data$ident)))) +
    ggplot2::coord_fixed() +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(size = 3))) +
    cowplot::theme_cowplot()
}

feature_scatter <- function(data, feature, reduction) {
  dim_names <- paste0(
    switch(reduction,
      tsne = "tSNE_", umap = "UMAP_", pca = "PC_"),
    1:2
  )

  col_name <- names(data)
  col_name[match(feature, col_name)] <- "expression"
  names(data) <- col_name

  data %>%
    dplyr::arrange(expression) %>%
    ggplot2::ggplot(
        mapping = ggplot2::aes_string(
          x = dim_names[1], y = dim_names[2], color = "expression")) +
      ggplot2::geom_point(size = .2) +
      ggplot2::scale_colour_gradient(low = "lightgrey", high = "blue") +
      ggplot2::ggtitle(feature) +
      ggplot2::coord_fixed() +
      cowplot::theme_cowplot() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.title = ggplot2::element_blank()
      )
}

feature_violin <- function(data, feature) {
  Seurat:::SingleExIPlot(
    type = "violin",
    data = data[, feature, drop = FALSE],
    idents = data$ident,
    adjust = 1, pt.size = 0
  ) +
  ggplot2::guides(fill = "none")
}

expression_projection_page <- function(
  input, output, session, library_list, cache) {
  shiny::updateSelectInput(
    session,
    "select_sample",
    choices = basename(names(library_list)),
  )

  gene_queries <- shiny::eventReactive(input$submit, {
    gene_list <- strsplit(stringr::str_trim(input$gene_list), "\n")[[1]]
    unique(gene_list)
  })

  # Loading Cell Embeddings
  cell_embeddings <- shiny::reactive({
    shiny::req(input$select_sample)
    fetch_cell_embeddings(
      library_list[[input$select_sample]]$Seurat_Disk, cache)
  })

  # Loading Gene Expressions
  gene_expressions <- shiny::eventReactive(input$submit, {
    shiny::req(input$select_sample)
    library <- input$select_sample

    fetch_gene_expression(
      library_list[[library]]$Seurat_Disk,
      gene_queries(),
      cache
    )
  })

  data_to_plot <- shiny::reactive({
    reduction <- input$reduction
    dim_names <- paste0(
      switch(reduction,
        tsne = "tSNE_", umap = "UMAP_", pca = "PC_"),
      1:2
    )
    cbind(
      cell_embeddings()[, c(dim_names, "ident")],
      gene_expressions()
    )
  })

  output$dimplot <- shiny::renderPlot({
    shiny::req(cell_embeddings())
    dim_names <- paste0(
      switch(input$reduction,
        tsne = "tSNE_", umap = "UMAP_", pca = "PC_"),
      1:2
    )

    data_to_plot <- cell_embeddings()[, c(dim_names, "ident")]
    cluster_plot(data_to_plot, reduction = input$reduction)
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

    invisible(lapply(plot_output_id_list, function(gene_id) {
      output[[paste0("scatter-", gene_id)]] <- shiny::renderPlot({
        data_to_plot <- data_to_plot()
        gene_id_p <- paste0("rna_", gene_id)
        if (gene_id_p %in% names(data_to_plot)) {
          feature_scatter(data_to_plot, gene_id_p, reduction = reduction)
        } else {
          plot_placeholder("Not Found")
        }
      })
      output[[paste0("violin-", gene_id)]] <- shiny::renderPlot({
        gene_id_p <- paste0("rna_", gene_id)
        data_to_plot <- data_to_plot()
        if (gene_id_p %in% names(data_to_plot)) {
          feature_violin(data_to_plot, gene_id_p)
        } else {
          plot_placeholder("Not Found")
        }
      })
    }))
  })
}