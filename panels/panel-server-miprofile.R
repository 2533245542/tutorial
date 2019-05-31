miprofile_UI_exist = FALSE
# this function checks the mutual exclusiveness of confouding and strata. If some strata other than "NULL" is chosen, 
# make confouding to have nothing selectd. Similarly, if some confouding is/are selected, make strat to have "NULL" selected
observeEvent(c(input$miprofile_conf), { 
    updateSelectInput(session, "miprofile_strata",selected = "NULL")
})
observeEvent(c(input$miprofile_strata), { 
    if(input$miprofile_strata != "NULL"){
      updateSelectInput(session, "miprofile_conf",selected = "")  
    }    
})
observeEvent(c(input$miprofile_uwu, input$miprofile_gu, input$miprofile_pwu), { # observe all three checkbox inputs
    
    if(!miprofile_UI_exist && sum(input$miprofile_uwu, input$miprofile_gu, input$miprofile_pwu) == 1){ # if no UI exists and 
      #only one checkboxinput is checked
      insertUI(
        # insert inside the div input
        selector = "#miprofile_fileinput_selector",
        where = "afterBegin",
        ui = fileInput("miprofile_treefile", "Upload .tre file", multiple = FALSE)
      )  
      miprofile_UI_exist <<- TRUE
    }
    
    if(miprofile_UI_exist && all(!input$miprofile_uwu, !input$miprofile_gu, !input$miprofile_pwu)){ # if UI exists and 
      # all are false, i.e. no one is checked, then remove UI    
      # remove the div that is inside the div with id fileinput we created
      removeUI(selector = "div#miprofile_fileinput_selector > div")
      miprofile_UI_exist <<- FALSE
    }
    
  })
# reactive functions
mistudio_miprofile_cova_result = reactive({

  source("functions/mistudio_milineage_cova.R")

  mistudio_milineage_cova(physeq())

})

mistudio_miprofile_conf_result = reactive({

  source("functions/mistudio_milineage_conf.R")

  mistudio_milineage_conf(physeq())

})

# reactive .R functions
mistudio_miprofile_result = reactive({

  source("functions/mistudio_miprofile.R")

  if(input$action_miprofile == 0){ # return NULL if run button unclicked
    return(NULL)
  }

  isolate({
    if(input$miprofile_gu) {
      gUniFrac.alpha = scan(text = input$miprofile_gu_alpha, what  = numeric(), sep = ",")
      gUniFrac.rarefy = rep(as.integer(input$miprofile_gu_norm), length(gUniFrac.alpha))
    } else {
      gUniFrac.alpha = NULL
    }

    if(input$miprofile_pwu) {      
      pUniFrac.alpha = scan(text = input$miprofile_pwu_alpha, what  = numeric(), sep = ",")
      pUniFrac.rarefy = rep(as.integer(input$miprofile_pwu_norm), length(pUniFrac.alpha))
    } else {
      pUniFrac.alpha = NULL
    }

    if(input$miprofile_uwu) {
      uwUniFrac.rarefy = as.integer(input$miprofile_uwu_norm)
    } else {
      uwUniFrac.rarefy = NULL
    }

    if(input$miprofile_bcd) {
      BC.rarefy = as.integer(input$miprofile_bcd_norm)
    } else {
      BC.rarefy = NULL
    }

    if(input$miprofile_jcd) {
      Jaccard.rarefy = as.integer(input$miprofile_jcd_norm)
    } else {
      Jaccard.rarefy = NULL
    }

    if(input$miprofile_strata == "NULL"){
      strata = NULL
    } else {
      strata = input$miprofile_strata
    }
    
    tree.file = input$miprofile_treefile$datapath

    result = mistudio_miprofile(
      phyobj= physeq(), 
      cova= input$miprofile_cova, conf = input$miprofile_conf, 
      gUniFrac.alpha = gUniFrac.alpha, gUniFrac.rarefy = gUniFrac.rarefy, 
      pUniFrac.alpha = pUniFrac.alpha, pUniFrac.rarefy = pUniFrac.rarefy, 
      uwUniFrac.rarefy = uwUniFrac.rarefy,
      BC.rarefy = BC.rarefy, 
      Jaccard.rarefy = Jaccard.rarefy,
      tree.file = tree.file,
      strata = strata,
      nperm = input$miprofile_nperm)

    return(result)
  })
})

mistudio_miprofile_pca_result = reactive({

  source("functions/mistudio_miprofile_pca.R")

  if(input$action_miprofile_pca == 0){ # return NULL if run button unclicked
    return(NULL)
  }

  isolate({
    if(input$miprofile_gu) {
      gUniFrac.alpha = scan(text = input$miprofile_gu_alpha, what  = numeric(), sep = ",")
      gUniFrac.rarefy = rep(as.integer(input$miprofile_gu_norm), length(gUniFrac.alpha))
    } else {
      gUniFrac.alpha = NULL
    }

    if(input$miprofile_pwu) {      
      pUniFrac.alpha = scan(text = input$miprofile_pwu_alpha, what  = numeric(), sep = ",")
      pUniFrac.rarefy = rep(as.integer(input$miprofile_pwu_norm), length(pUniFrac.alpha))
    } else {
      pUniFrac.alpha = NULL
    }

    if(input$miprofile_uwu) {
      uwUniFrac.rarefy = as.integer(input$miprofile_uwu_norm)
    } else {
      uwUniFrac.rarefy = NULL
    }

    if(input$miprofile_bcd) {
      BC.rarefy = as.integer(input$miprofile_bcd_norm)
    } else {
      BC.rarefy = NULL
    }

    if(input$miprofile_jcd) {
      Jaccard.rarefy = as.integer(input$miprofile_jcd_norm)
    } else {
      Jaccard.rarefy = NULL
    }

    if(input$miprofile_color == "NULL"){
      color = NULL
    } else {
      color = input$miprofile_color
    }
    if(input$miprofile_shape == "NULL"){
      shape = NULL
    } else {
      shape = input$miprofile_shape
    }
    # if(input$miprofile_facetrow == "NULL"){
    #   facetrow = NULL
    # } else {
    #   facetrow = input$miprofile_facetrow
    # }
    # if(input$miprofile_facetcol == "NULL"){
    #   facetcol = NULL
    # } else {
    #   facetcol = input$miprofile_facetcol
    # }
    if(input$miprofile_label == "NULL"){
      label = NULL
    } else {
      label = input$miprofile_label
    }
    
    tree.file = input$miprofile_treefile$datapath

    result = mistudio_miprofile_pca(
      phyobj= physeq(), 
      gUniFrac.alpha = gUniFrac.alpha, gUniFrac.rarefy = gUniFrac.rarefy, 
      pUniFrac.alpha = pUniFrac.alpha, pUniFrac.rarefy = pUniFrac.rarefy, 
      uwUniFrac.rarefy = uwUniFrac.rarefy,
      BC.rarefy = BC.rarefy, 
      Jaccard.rarefy = Jaccard.rarefy,
      tree.file = tree.file,
      nperm = input$miprofile_nperm,
      color = color,
      shape = shape,
      facetrow = input$miprofile_facetrow,
      facetcol = input$miprofile_facetcol,
      label = label,
      label_size = input$miprofile_label_size,
      label_x = input$miprofile_label_x,
      label_y = input$miprofile_label_y,
      label_rotate = input$miprofile_label_rotate,      
      palette = input$miprofile_palette,
      theme = input$miprofile_theme,
      point_size = input$miprofile_point_size,
      point_opacity = input$miprofile_point_opacity,
      xlab = input$miprofile_xlab,
      ylab = input$miprofile_ylab,
      width = input$miprofile_width,
      height = input$miprofile_height
      )
    return(result)
  })
})

# UI outputs
output$miprofile_out_cova <- renderUI({
  req(physeq)
  selectInput("miprofile_cova", "Covariate", choices = mistudio_miprofile_cova_result(), multiple = TRUE) # although UI is changed, the server side is not supported for multiple cova processing
})

output$miprofile_out_conf <- renderUI({
  req(physeq)
  selectInput("miprofile_conf", "Confouding", choices = mistudio_miprofile_conf_result(), multiple = TRUE)
})

output$miprofile_out_strata <- renderUI({
  req(physeq)
  selectInput("miprofile_strata", "Strata", choices = c("NULL",sample_variables(physeq())), multiple = FALSE)
})


output$miprofile_out_color <- renderUI({
  selectInput("miprofile_color", "Color", choices = c("NULL",sample_variables(physeq())))
})
output$miprofile_out_shape <- renderUI({
  selectInput("miprofile_shape", "Shape", choices = c("NULL",sample_variables(physeq())))
})
output$miprofile_out_facetrow <- renderUI({
  selectInput("miprofile_facetrow", "Facet Row", choices = sample_variables(physeq()), multiple = TRUE)
})
output$miprofile_out_facetcol <- renderUI({
  selectInput("miprofile_facetcol", "Facet Col", choices = sample_variables(physeq()), multiple = TRUE)
})
output$miprofile_out_label <- renderUI({
  selectInput("miprofile_label", "Label", choices = c("NULL",sample_variables(physeq())))
})



# outputs
output$miprofile_pval = renderTable({
  req(mistudio_miprofile_result())
  result = mistudio_miprofile_result()
  result = as.data.frame(t(result))
  return(result)
}, digits = 10)
output$miprofile_pca_ui = renderUI({  
    plotOutput("miprofile_pca", width = input$miprofile_width, height = input$miprofile_height)
})
output$miprofile_pca = renderPlot({
    req(mistudio_miprofile_pca_result())
    result = mistudio_miprofile_pca_result()
    return(result)
})
    


outputOptions(output, "miprofile_pca", suspendWhenHidden = FALSE)  
