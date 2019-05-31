library(ggplot2)
source("panels/panel-ui-file.R", local = TRUE)
source("panels/panel-ui-filter.R", local = TRUE)
source("panels/panel-ui-milineage.R", local = TRUE)
source("panels/panel-ui-miprofile.R", local = TRUE)
source("panels/panel-ui-milineagec.R", local = TRUE)

ui = navbarPage(
  title = h4(a(href="https://tangzheng1.github.io/tanglab/software.html",   "miStudio")),
  tabPanel("File", file),
  tabPanel("Filter", filter),
  tabPanel("miLineage", milineage),
  tabPanel("milineageC", milineagec),
  tabPanel("miProfile", miprofile),
  # inverse = TRUE,
  collapsible = TRUE,
  header = headerTagList,
  windowTitle = "miStudio"
)
shinyUI(ui)