################################################################################
shiny_busy <- function() {
  # use &nbsp; for some alignment, if needed
  HTML("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;", paste0(
    '<span data-display-if="',
    '$(&#39;html&#39;).attr(&#39;class&#39;)==&#39;shiny-busy&#39;',
    '">',
    '<i class="fa fa-spinner fa-pulse fa-fw" style="color:orange"></i>',
    '</span>'
  ))
}

miprofile = fluidPage(
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(type = "tabs",
                        tabPanel("Analysis",
                          actionButton("action_miprofile", "RUN", icon("filter")),
                          uiOutput("miprofile_out_cova"),
                          uiOutput("miprofile_out_conf"),
                          uiOutput("miprofile_out_strata"),

                          fluidRow(style = "border: 2px double gray;",
                            column(width = 5, div(style="font-size:125%;color:black;",  h5("Distances"))),
                            column(width = 4, div(style="font-size:125%;color:black;",  h5("Alpha(s)"))),
                            column(width = 3, div(style="font-size:125%;color:black;",  h5("Rarefy?")))        
                          ),
                          fluidRow(
                            # column(width = 6, checkboxInput("miprofile_bcd", "Bray-Curtis Distance", value = TRUE)),
                            column(width = 9, div(style="font-size:100%;color:black;",checkboxInput("miprofile_bcd", "Bray-Curtis", value = TRUE))),
                            
                            column(width = 3, checkboxInput("miprofile_bcd_norm", label = NULL, value = FALSE, width = NULL)) # if true, means rarefy, as.integer(TRUE) returns 1, so wouldn't affect the server side       

                            # column(width = 3, selectInput("miprofile_bcd_norm", label = NULL, choices = list(Y = 1, N = 0), selected = 0))        
                          ),

                          # div(style="font-size:125%;color:blue;", checkboxInput("miprofile_bcd", "Bray-Curtis Distance", value = TRUE)),      
                          # selectInput("miprofile_bcd_norm", "Rarefy?",list(Yes = 1, No = 0), selected = 0),
                          
                          fluidRow(
                            column(width = 9, div(style="font-size:100%;color:black;",checkboxInput("miprofile_jcd", "Jaccard", value = TRUE))),
                            column(width = 3, checkboxInput("miprofile_jcd_norm", label = NULL, value = TRUE, width= NULL))        
                          ),

                          # div(style="font-size:125%;color:blue;",checkboxInput("miprofile_jcd", "Jaccard Distance", value = TRUE)),
                          # selectInput("miprofile_jcd_norm", "Rarefy?",list(Yes = 1, No = 0), selected = 1),

                          div(id = "miprofile_fileinput_selector"), # for displaying file input handler dynamically

                          # fileInput("miprofile_treefile", "Upload .tre file", multiple = FALSE),
                          
                          fluidRow(
                            column(width = 4, div(style="font-size:100%;color:black;",checkboxInput("miprofile_uwu", "UnWeighted UniFrac", value = FALSE))),
                            column(width = 5),
                            column(width = 3, checkboxInput("miprofile_uwu_norm", label = NULL , value = TRUE, width= NULL))        
                          ),

                          # div(style="font-size:125%;color:blue;",checkboxInput("miprofile_uwu", "UnWeighted UniFrac", value = FALSE)),
                          # selectInput("miprofile_uwu_norm", "Rarefy?",list(Yes = 1, No = 0), selected = 0),
                          fluidRow(
                            column(width = 4, div(style="font-size:100%;color:black;",checkboxInput("miprofile_gu", "Generalized UnWeighted UniFrac", value = FALSE))),
                            column(width = 5, textInput("miprofile_gu_alpha", label = NULL, value = "0")),
                            column(width = 3, checkboxInput("miprofile_gu_norm", label = NULL , value = FALSE, width= NULL))        
                          ),

                          fluidRow(
                            column(width = 4, div(style="font-size:100%;color:black;",checkboxInput("miprofile_pwu", "Presence Weighted UniFrac", value = FALSE))),
                            column(width = 5, textInput("miprofile_pwu_alpha", label = NULL, value = "0")),
                            column(width = 3, checkboxInput("miprofile_pwu_norm", label = NULL , value = TRUE, width= NULL))        
                          ),
                          # div(style="font-size:125%;color:blue;",checkboxInput("miprofile_gu", "Generalized UniFrac", value = FALSE)),      
                          # textInput("text1", label = "Dist Norm", value = "P"),
                          # textInput("text2", label = "Alpha", value = "0"),
                          
                          # div(style="font-size:125%;color:blue;",checkboxInput("miprofile_pwu", "Presence Weighted UniFrac", value = FALSE)),
                          # textInput("text3", label = "Dist Norm", value = "P"),
                          # textInput("text4", label = "Alpha", value = "0"),

                          div(style="font-size:90%;color:gray;",h5("Alpha(s) support multi-value input seperated by comma, e.g. Alpha(s) = 0,0.1,0.2,0.3")),
                          helpText("Alpha(s) support multi-value input seperated by comma, e.g. Alpha(s) = 0,0.1,0.2,0.3"),
                          numericInput("miprofile_nperm", "nperm", 1000000, step = 100000)
                        ),

                        tabPanel("Visualization",
                          actionButton("action_miprofile_pca", "Plot", icon("signal")),
                          fluidRow(
                            column(width = 6, uiOutput("miprofile_out_color")),
                            column(width = 6, uiOutput("miprofile_out_shape"))
                          ),
                          fluidRow(
                            column(width = 6, uiOutput("miprofile_out_facetrow")),
                            column(width = 6, uiOutput("miprofile_out_facetcol"))
                          ),
                          fluidRow(
                            column(width = 6, uiOutput("miprofile_out_label")),
                            column(width = 6, numericInput("miprofile_label_size", "Label size", value = 3, step = 1))
                          ),
                          fluidRow(                            
                            column(width = 4, numericInput("miprofile_label_x", "Move lable up/down by", value = 0, step = 1)),
                            column(width = 4, numericInput("miprofile_label_y", "Move lable left/right by", value = 0, step = 1)),
                            column(width = 4, numericInput("miprofile_label_rotate", "Rotate lable by ", value = 0, step = 10))
                          ),

                          fluidRow(
                            column(width = 6, selectInput("miprofile_palette", "Select palette", choices = list(default = "Set1", set2 = "Set2", set3 = "Set3", blue = 1, green = 2, purple = 3, gray =  6, orange = 7, pink = 13))),                            
                            column(width = 6, selectInput("miprofile_theme", "Select theme", choices = list(default = "theme_bw", thin = "theme_linedraw", light = "theme_light",minimal = "theme_minimal",classic = "theme_classic",gray = "theme_gray")))                            
                          ),

                          fluidRow(
                            column(width = 6, numericInput("miprofile_point_size", "Point size", value = 1, step = 1)),                            
                            column(width = 6, numericInput("miprofile_point_opacity", "Point opacity", value = 1, step = 0.1, min = 0, max = 1))
                          ),

                          fluidRow(
                            column(width = 6, textInput("miprofile_xlab", "X lab")),                            
                            column(width = 6, textInput("miprofile_ylab", "Y lab"))
                          ),                          

                          fluidRow(
                            column(width = 4, numericInput("miprofile_width", "Width", value = 500, step = 50)),
                            column(width = 4, numericInput("miprofile_height", "Height", value = 500, step = 50)),                            
                            column(width = 4, selectInput("miprofile_format", "Format", choices = c("pdf", "png", "jpeg"), multiple = FALSE))
                          )
                                                                  
                        ),
                        tabPanel("Download",
                          downloadButton("miprofile_download_work", "Downlaod Work")
                          )


      )
      
    ),

    # Main Panel
    mainPanel(
      tabsetPanel(type = "pills",
                        tabPanel("p-value",
                          shiny_busy(),                                          
                          htmlOutput("miprofile_pval")
                          ),
                        tabPanel("PCA plots", uiOutput("miprofile_pca_ui"))  
      )
    )
    # column(
    #   width = 8, offset = 0, 
    #   h4("miProfile Result"),
    #   h5("p-value"),
    #   shiny_busy(),
      
    #   plotOutput("miprofile_pca_plot")
    # )
  )
)
################################################################################
