options(shiny.maxRequestSize = 100*1024^2)
# Set Shiny Reaction Log to TRUE
options(shiny.reactlog=TRUE)
# Default ggplot2 theme (Only relevant if panel-specific theme missing or NULL)
# theme_set(theme_bw())
# library(ggplot2)
shinyServer(function(input, output){
  source("panels/panel-server-file.R", local = TRUE)
  source("panels/panel-server-filter.R", local = TRUE)
  source("panels/panel-server-milineage.R", local = TRUE)
  source("panels/panel-server-miprofile.R", local = TRUE)

})