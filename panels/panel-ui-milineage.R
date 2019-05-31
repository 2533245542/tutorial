################################################################################

# library(shiny)

milineage = fluidPage(
  sidebarLayout(
    # Sidebar Panel
    sidebarPanel(
      tabsetPanel(type = "tabs",
                        tabPanel("Analysis",
                          selectInput("milineage_function", "Select function to run", choices = c("QCAT", "QCAT_GEE", "ZIGDM")),
                          uiOutput("milineage_out_cova"),
                          div(id = "milineage_cova_selector"), # for displaying QCAT_GEE select dynamically
                          uiOutput("milineage_out_conf"),
                          div(id = "milineage_conf_selector"), # for displaying QCAT_GEE and ZIGDM select dynamically                                                    

                          fluidRow(
                            column(width = 4, numericInput("milineage_mindepth", "Min depth", 0)),
                            column(width = 4, numericInput("milineage_nresample", "n.resample", 1)),                            
                            column(width = 4, numericInput("milineage_fdralpha", "Fdr.alpha", 0.05))
                          ),
                          
                          fluidRow(
                            column(width = 6, div(id = "milineage_ZILB_selector")),  # for displaying ZIGDM select dynamically 
                                                                                     # CANNOT use ZI.LB for selector cuz syntax conflict                         
                            column(width = 6, div(id = "milineage_testtype_selector"))  # for displaying ZIGDM select dynamically                          
                          ),
                          actionButton("action_milineage", "Run", icon("filter"))
                        ),

                        tabPanel("Visualization", 
                          # uiOutput("milineage_out_lineage"), 
                          uiOutput("milineage_out_annot_lineage"), 
                          uiOutput("milineage_out_compositional_lineage"), 
                          selectInput("milineage_arrange", "Type of x axis", choices = c(Categorical = "categorical", continuous = "continuous"), multiple = FALSE),                        
                          fluidRow(
                            column(width = 6, div(id = "milineage_categorical_selector")),
                            column(width = 6, div(id = "milineage_subcategorical_selector"))
                          ),
                          div(id = "milineage_continuous_selector"),

                          fluidRow(
                            column(width = 6, uiOutput("milineage_out_stratify")),
                            column(width = 6, uiOutput("milineage_out_substratify"))                                                    
                          ),

                          fluidRow(
                            column(width = 6, selectInput("milineage_palette", "Palette", choices = list(default = "default", set1 = "Set1", set2 = "Set2", set3 = "Set3", blue = 1, green = 2, purple = 3, gray =  6, orange = 7, pink = 13))),                            
                            column(width = 6, selectInput("milineage_theme", "Theme", choices = list(default = "default", thin = "theme_linedraw", light = "theme_light",minimal = "theme_minimal",classic = "theme_classic",gray = "theme_gray")))                            
                          ),
                          
                          fluidRow(
                            column(width = 4, textInput("milineage_xlab", "X lab")),                            
                            column(width = 4, textInput("milineage_ylab", "Y lab")),
                            column(width = 4, textInput("milineage_title", "Title"))                            
                          ),
                          fluidRow(                          
                            column(width = 6, numericInput("milineage_proportion", "Minimimum lineage proportion", value = 0.01, min = 0, max = 1, step = 0.01)),
                            column(width = 6, checkboxInput("milineage_hidex", "Hide X labels and ticks"))                            
                          ),
                          #fluidRow(
                          #  column(width = 6, uiOutput("milineage_out_sortby")),
                          #  column(width = 6, uiOutput("milineage_out_stratifyby"))
                          #),
                          fluidRow(
                            column(width = 4, numericInput("milineage_width", "Width", value = 500, step = 50)),
                            column(width = 4, numericInput("milineage_height", "Height", value = 500, step = 50)),                            
                            column(width = 4, selectInput("milineage_format", "Format", choices = c("pdf", "png", "jpeg"), multiple = FALSE))
                          ),                                              
                          actionButton("action_milineage_compositional", "Plot", icon = icon("signal"))

                        ),
                        tabPanel("Download", 
                          downloadButton("milineage_reproduce", "Downlaod Work")
                          )

      )
      
    ),
    mainPanel(
      # Main Panel
      tabsetPanel(type = "pills",
                        tabPanel("Lineage p-value", 
                          tableOutput("lineage.pval"),
                          tableOutput("global.pval"),
                          tableOutput("sig.lineage")
                          ),
                        tabPanel("Tree Annotation", plotOutput("milineage_tree")),
                        tabPanel("Compositional Plot", plotOutput("milineage_compositional"))
      )
    )
  )
)
################################################################################
