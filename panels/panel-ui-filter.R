################################################################################
filter = fluidPage(
  fluidRow(
    # Sidebar Panel default is 4-column.
    sidebarPanel(

      # selectInput("filter_subset_ranks", "Subset Ranks", list("NULL")),
      # selectInput("filter_subset_taxa", "Subset Taxa",list("NULL")),
      uiOutput("filter_uix_subset_taxa_ranks"),
      uiOutput("filter_uix_subset_taxa_select"),
      uiOutput("filter_uix_subset_sample_vars"),
      uiOutput("filter_uix_subset_sample_select"),

      # selectInput("filter_subset_variables", "Subset sample variables",list("NULL")),
      # selectInput("filter_subset_classes", "Subset variable classes",list("NULL")),
      # numericInput("filter_OTUsums", "Select sample minimum OTU sums", 1000),
      # numericInput("filter_taxasums", "Select taxa minimum counts", 10)
      numericInputRow("filter_sample_sums_threshold", "Select sample minimum OTU sums",
                                          value=SampleSumDefault, min=0, step=100, class="col-md-12"),
      numericInputRow("filter_taxa_sums_threshold", "Select taxa minimum counts",
                                          value=OTUSumDefault, min=0, step=1, class="col-md-12"),
      
      actionButton("action_filter", "Filter", icon("filter")),
      actionButton("action_export", "Export Data", icon("filter"))
    ),
    # Now the Main Panel.
    column(
      width = 8, offset = 0, 
      h4("Data Summaries"),
      fluidRow(
        column(width = 6,
               h4("Original"),
               htmlOutput('filtered_contents0')
        ),
        column(width = 6,
               h4("Filtered Data:"),
               htmlOutput('filtered_contents')               
      )),
      h4("Component Table, Filtered Data"),
      fluidRow(column(width=12,
          div(class="col-md-8", uiOutput("uix_available_components_filt")),
          div(class="col-md-3", numericInputRow("component_table_colmax_filt", "Max. Columns",
                                             value = 25L, min = 1L, step = 5L, class="col-md-12"))
      )),
      dataTableOutput('physeqComponentTable')
    )
  )
)
################################################################################
