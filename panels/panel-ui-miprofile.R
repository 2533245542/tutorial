################################################################################
miprofile = fluidPage(
  fluidRow(
    # Sidebar Panel default is 4-column.
    textOutput("test"),
    sidebarPanel(
      actionButton("action_miprofile", "RUN", icon("filter")),
      uiOutput("miprofile_out_cova"),
      uiOutput("miprofile_out_conf"),
      h4("Warning: Do not select identical variables for COVA and CONF"),
      
      div(style="font-size:125%;color:blue;",checkboxInput("miprofile_bcd", "Bray-Curtis Distance", value = TRUE)),
      selectInput("miprofile_bcd_norm", "DIST_NORM",list(Rarefaction = "Rarefaction", Proportion = "Proportion"), selected = "Proportion"),
      
      div(style="font-size:125%;color:blue;",checkboxInput("miprofile_jcd", "Jaccard Distance", value = TRUE)),
      selectInput("miprofile_jcd_norm", "DIST_NORM",list(Rarefaction = "Rarefaction", Proportion = "Proportion"), selected = "Rarefaction"),

      
      div(style="font-size:125%;color:blue;",checkboxInput("miprofile_uwu", "UnWeighted UniFrac", value = TRUE)),
      selectInput("miprofile_uwu_norm", "DIST_NORM",list(Rarefaction = "Rarefaction", Proportion = "Proportion"), selected = "Rarefaction"),
      
      div(style="font-size:125%;color:blue;",checkboxInput("miprofile_gu", "Generalized UniFrac", value = TRUE)),
      selectInput("miprofile_gu_norm", "DIST_NORM",list(Rarefaction = "Rarefaction", Proportion = "Proportion"), selected = "Proportion"),
      numericInput("miprofile_gu_alpha", "alpha", 0),
      
      div(style="font-size:125%;color:blue;",checkboxInput("miprofile_pwu", "Presence Weighted UniFrac", value = TRUE)),
      selectInput("miprofile_pwu_norm", "DIST_NORM",list(Rarefaction = "Rarefaction", Proportion = "Proportion"), selected = "Rarefaction"),
      numericInput("miprofile_pwu_alpha", "alpha", 0),

      fileInput("filedistance", "Upload DIST_IFILE", multiple = TRUE),

      checkboxInput("miprofile_outfiles", "OUTPUT DISTANCES INTO FILES", value = TRUE)

    )
  )
)
################################################################################
