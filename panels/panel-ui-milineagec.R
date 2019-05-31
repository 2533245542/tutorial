milineagec = fluidPage(
  sidebarPanel(
    tabsetPanel(type = "tabs",
                      tabPanel("Analysis",
                        selectInput("milineagec_function", "Select function to run", choices = c("QCAT.Cluster", "QCAT_GEE.Cluster")),
                        uiOutput("milineagec_out_cova"),
                        div(id = "milineagec_cova_selector"), # for displaying QCAT_GEE select dynamically
                        uiOutput("milineagec_out_conf"),
                        div(id = "milineagec_conf_selector"), # for displaying QCAT_GEE and ZIGDM select dynamically                          
                        fluidRow(
                          column(width = 6, uiOutput("milineagec_out_id")),
                          column(width = 6, uiOutput("milineagec_out_timebase"))                            
                        ),
                        fluidRow(
                          column(width = 6, uiOutput("milineagec_out_permtype")),
                          column(width = 6, uiOutput("milineagec_out_testtype"))                                                    
                        ), 
                        fluidRow(
                          column(width = 4, numericInput("milineagec_mindepth", "Min depth", 0)),
                          column(width = 4, numericInput("milineagec_nresample", "n.resample", 1000)),                            
                          column(width = 4, numericInput("milineagec_fdralpha", "Fdr.alpha", 0.05))
                        ),   
                        actionButton("action_milineagec", "Run", icon("filter"))
                      ),

                      tabPanel("Visualization",
                        actionButton("action_milineagec_compositional", "Plot compositinal", icon("filter")),

                        uiOutput("milineagec_out_compositional_lineage"),
                        selectInput("milineagec_arrange", "Type of x axis", choices = c(Categorical = "categorical", continuous = "continuous"), multiple = FALSE),                        
                        fluidRow(
                          column(width = 6, div(id = "milineagec_categorical_selector")),
                          column(width = 6, div(id = "milineagec_subcategorical_selector"))
                        ),
                        div(id = "milineagec_continuous_selector"),

                        fluidRow(
                          column(width = 6, uiOutput("milineagec_out_stratify")),
                          column(width = 6, uiOutput("milineagec_out_substratify"))                                                    
                        ),

                        fluidRow(
                          column(width = 6, selectInput("milineagec_palette", "Pallete", choices = list(default = "default", set1 = "Set1", set2 = "Set2", set3 = "Set3", blue = 1, green = 2, purple = 3, gray =  6, orange = 7, pink = 13))),                            
                          column(width = 6, selectInput("milineagec_theme", "Theme", choices = list(default = "default", thin = "theme_linedraw", light = "theme_light",minimal = "theme_minimal",classic = "theme_classic",gray = "theme_gray")))                            
                        ),
                    
                        fluidRow(
                          column(width = 4, textInput("milineagec_xlab", "X lab")),                            
                          column(width = 4, textInput("milineagec_ylab", "Y lab")),
                          column(width = 4, textInput("milineagec_title", "Title"))                            
                        ),
                        fluidRow(                          
                          column(width = 6, numericInput("milineagec_proportion", "Minimimum lineage proportion", value = 0.01, min = 0, max = 1, step = 0.01)),
                          column(width = 6, checkboxInput("milineagec_hidex", "Hide X lables and ticks"))                            
                        ),
                        fluidRow(
                          column(width = 4, numericInput("milineagec_width", "Width", value = 500, step = 50)),
                          column(width = 4, numericInput("milineagec_height", "Height", value = 500, step = 50)),                            
                          column(width = 4, selectInput("milineagec_format", "Format", choices = c("pdf", "png", "jpeg"), multiple = FALSE))
                        )                                              
                      ),
                      tabPanel("Download",
                        downloadButton("milineagec_compositional_reproduce", "Download Work")
                        )

    )
    
  ),
  
  mainPanel(
  	tabsetPanel(type = "pills",
  	                      tabPanel("Result", 
                            tableOutput("milineagec_lineagepval"),
                            tableOutput("milineagec_globalpval"),
                            tableOutput("milineagec_siglineage")                    
  	                      ),

  	                      tabPanel("plot",
                            plotOutput("milineagec_compositional")                       
  	                      )
  	)
  	# tableOutput("filetree_info")
  )
)
