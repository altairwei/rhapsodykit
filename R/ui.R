#' Create main UI for application
#'
create_main_ui <- function() {
  shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(
      title = "Rhapsody WTA Viewer"
    ),
    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        shinydashboard::menuItem("Metrics Summary", tabName = "metrics_summary")
      )
    ),
    shinydashboard::dashboardBody(
      shinydashboard::tabItems(
        shinydashboard::tabItem(
          tabName = "metrics_summary",
          metrics_summary_page_ui("metrics_summary_page")
        )
      )
    )
  )
}