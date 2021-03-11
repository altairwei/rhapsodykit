#' Create main UI for application
#'
create_main_ui <- function() {
  shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(
      title = "Rhapsody WTA Viewer"
    ),
    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        # shinydashboard::menuItem(
        #   "Metrics Summary", tabName = "metrics_summary"),
        shinydashboard::menuItem("Standard Analysis",
          shinydashboard::menuSubItem(
            "Expression Projection", tabName = "expression_projection"),
          tabName = "standard_analysis"
        ),
        shinydashboard::menuItem("Utilities",
          shinydashboard::menuSubItem(
            "Ensembl Homologs", tabName = "ensembl_homologs"),
          tabName = "utilities"
        )
      )
    ),
    shinydashboard::dashboardBody(
      waiter::use_waiter(),
      waiter::use_waitress(),
      waiter::waiter_show_on_load(html = waiter::spin_flower()),
      shinydashboard::tabItems(
        # shinydashboard::tabItem(
        #   tabName = "metrics_summary",
        #   metrics_summary_page_ui("metrics_summary_page")
        # ),
        shinydashboard::tabItem(
          tabName = "expression_projection",
          expression_projection_page_ui("expression_projection_page")
        ),
        shinydashboard::tabItem(
          tabName = "ensembl_homologs",
          ensembl_homologs_page_ui("ensembl_homologs_page")
        )
      )
    )
  )
}