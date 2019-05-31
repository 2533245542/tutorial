################################################################################
filter = fluidPage(
  fluidRow(
    # Sidebar Panel default is 4-column.
    sidebarPanel(

      # selectInput("filter_subset_ranks", "Subset Ranks", list("NULL")),
      # selectInput("filter_subset_taxa", "Subset Taxa",list("NULL")),
      uiOutput("filter_out_otu_aggregate"),
      fluidRow(
        column(width = 6, uiOutput("filter_out_subset_rank")),
        column(width = 6, uiOutput("filter_out_subset_taxa"))
        ),
      fluidRow(
        column(width = 6, uiOutput("filter_out_subset_subject")),
        column(width = 6, uiOutput("filter_out_subset_variable"))
        ),
      # fluidRow(
      #   column(width = 6, uiOutput("filter_uix_subset_taxa_ranks")),
      #   column(width = 6, uiOutput("filter_uix_subset_taxa_select"))
      #   ),
      
      # fluidRow(
      #   column(width = 6, uiOutput("filter_uix_subset_sample_vars")),
      #   column(width = 6, uiOutput("filter_uix_subset_sample_select"))
      #   ),
      fluidRow(
        column(width = 6, numericInputRow("filter_otu_minsum", "Minimum OTU count",
                                          value=0, min=0, step=5000, class="col-md-12")),
        column(width = 6, numericInputRow("filter_otu_zerofraction", "Maximum zero proportion",
                                          value=1.0, min=0, max = 1, step=0.1, class="col-md-12"))
        ),      
      actionButton("action_filter", "Filter", icon("filter")),
      downloadButton("action_export", "Export Data")
      # actionButton("action_export", "Export Data", icon("filter"))
    ),
    # Now the Main Panel.
    column(
      width = 8, offset = 0, 
      h4("Data Summaries"),
      fluidRow(
        h4("Before filtering"),
        column(width = 6,
               h4("hist(rowSums(otu))"),
               plotOutput('filter_plot_min_sum_before_filter')
        ),
        column(width = 6,
               h4("hist(the proportion of zero in each otu)"),

               plotOutput('filter_plot_zerofraction_before_filter')               
      )),
      fluidRow(
        h4("After filtering"),
        column(width = 6,
               h4("hist(rowSums(otu_filter))"),
               plotOutput('filter_plot_min_sum_after_filter')
        ),
        column(width = 6,
               h4("hist(the proportion of zero in each otu_filter)"),

               plotOutput('filter_plot_zerofraction_after_filter')               
      )),
      fluidRow(
        column(width = 6,
               h4("Original"),
               htmlOutput('filtered_contents0')
        ),
        column(width = 6,
               h4("Filtered Data:"),
               htmlOutput('filtered_contents')               
      )),
      h4("Tax table summary, original"),
      dataTableOutput("filter_summary_tax_before_filter"),
      h4("Tax table summary, filtered"),
      dataTableOutput("filter_summary_tax_after_filter"),
      h4("Meta data summary, original"),
      dataTableOutput("filter_summary_sample_before_filter"),
      h4("Meta data summary, filtered"),
      dataTableOutput("filter_summary_sample_after_filter"),
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
