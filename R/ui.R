#' Create main UI for application
#'
create_main_ui <- function() {
  shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(
      title = "Rhapsody WTA Viewer"
    ),
    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        shinydashboard::menuItem(
          "Metrics Summary", tabName = "metrics_summary"),
        shinydashboard::menuItem("Standard Analysis",
          shinydashboard::menuSubItem(
            "Expression Projection", tabName = "expression_projection"),
          tabName = "standard_analysis"
        )
      )
    ),
    shinydashboard::dashboardBody(
      waiter::use_waiter(),
      shinydashboard::tabItems(
        shinydashboard::tabItem(
          tabName = "metrics_summary",
          metrics_summary_page_ui("metrics_summary_page")
        ),
        shinydashboard::tabItem(
          tabName = "expression_projection",
          expression_projection_page_ui("expression_projection_page")
        )
      )
    )
  )
}