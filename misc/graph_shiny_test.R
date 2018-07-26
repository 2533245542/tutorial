library(shiny)

ui <- pageWithSidebar(
  headerPanel("renderImage example"),
  sidebarPanel(
    sliderInput("obs", "Number of observations:",
                min = 0, max = 1000,  value = 500)
  ),
  mainPanel(
    # Use imageOutput to place the image on the page
    imageOutput("myImage")
  )
)

server <- function(input, output, session) {
  output$myImage <- renderImage({
      filename = "migraph.png"
      list(src = filename)
  }, deleteFile = TRUE)
}

shinyApp(ui, server)