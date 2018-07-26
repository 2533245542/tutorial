################################################################################

milineage = fluidPage(
  fluidRow(
    # Sidebar Panel
    sidebarPanel(
      h4("miLineage"),
      fluidRow(
        column(width = 6, uiOutput("milineage_out_cov_n")),
        column(width = 6, uiOutput("milineage_out_cov_c"))
      ),
      fluidRow(
        column(width = 6, uiOutput("milineage_out_con_n")),
        column(width = 6, uiOutput("milineage_out_con_c"))
      ),

      fluidRow(
        column(width = 4, numericInput("milineage_mindepth", "Min depth", 0)),
        column(width = 4, numericInput("milineage_nresample", "n.resample", 1000)),
        column(width = 4, numericInput("milineage_fdralpha", "Fdr.alpha", 0.05))
      ),
      actionButton("action_milineage", "Run", icon("filter")),

      h4("Compositional Plot"),

      uiOutput("milineage_out_lineage"), 
      fluidRow(
        column(width = 6, uiOutput("milineage_out_sortby")),
        column(width = 6, uiOutput("milineage_out_stratifyby"))
      ),

      actionButton("action_milineage_compositional", "Plot", icon = icon("signal")),

      downloadButton("downloadData", "Downlaod Work")
    ),
    # Main Panel
    column(
      width = 8, offset = 0, 
      h4("QCAT_GEE Result"),
      h5("lineage p-value"),
      htmlOutput("lineage.pval"),
      h5("Global p-value"),
      htmlOutput("global.pval"),
      h5("Significant Lineage"),
      tableOutput("sig.lineage"),
      h4("Graph Representation"),
      imageOutput("myImage"),
      h4("Compositional Plot"),
      plotOutput("milineage_compositional")
    )
  )
)
################################################################################
