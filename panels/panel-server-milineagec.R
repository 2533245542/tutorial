
# make dynamically slection menu based on functions being selected
observeEvent({input$milineagec_function}, { # observe what function user selects
      
      if(input$milineagec_function == "QCAT_GEE.Cluster"){ 
        updateSelectInput(session, "milineagec_cova",
              label = "Covariate for positive part"
              )

        updateSelectInput(session, "milineagec_conf",
              label = "Confouding for positive part"
              )
        
        insertUI(
          selector = "#milineagec_cova_selector",
          where = "afterBegin",
          ui =  selectInput("milineagec_cova_extra", "Covariate for zero part", mistudio_milineagec_cova_result(), 
            selected = input$milineagec_cova, multiple = TRUE)
        )  
        insertUI(
          selector = "#milineagec_conf_selector",
          where = "afterBegin",
          ui =  selectInput("milineagec_conf_extra", "Confounding for zero part", mistudio_milineagec_conf_result(), 
            selected = input$milineagec_conf, multiple = TRUE)
        )    
        
      } else if(input$milineagec_function == "QCAT.Cluster"){
        updateSelectInput(session, "milineagec_cova",
              label = "Covariate"
              )

        updateSelectInput(session, "milineagec_conf",
              label = "Confouding"
              )
        removeUI(selector = "div#milineagec_cova_selector > div")
        removeUI(selector = "div#milineagec_conf_selector > div")

      } else{

      }

})

# let cova(conf) and cova_extra(conf_extra) always have the same value unless explicily modified
observeEvent(c(input$milineagec_cova, input$milineagec_conf), { # using c() is the correct way to observe two events

  updateSelectInput(session, "milineagec_cova_extra",
        selected = input$milineagec_cova
        )

  updateSelectInput(session, "milineagec_conf_extra",
        selected = input$milineagec_conf
        )
})

observeEvent({input$milineagec_arrange}, { # observe what function user selects
  if(is.null(physeq())){
    sample_col_choices = NULL
  } else {
    sample_col_choices = c("NULL", sample_variables(physeq()))
  }

  if(input$milineagec_arrange == "categorical"){ 

    removeUI(selector = "div#milineagec_continuous_selector > div")
    
  
    insertUI(
      selector = "#milineagec_categorical_selector",
      where = "afterBegin",
      ui =  selectInput("milineagec_categorical", "Group and average sample by", choices = sample_col_choices, 
        selected = NULL, multiple = FALSE)
    ) 
    insertUI(
      selector = "#milineagec_subcategorical_selector",
      where = "afterBegin",
      ui =  selectInput("milineagec_subcategorical", "Selected categories", choices = NULL, 
        selected = NULL, multiple = TRUE)
    )  

    
  } else if (input$milineagec_arrange == "continuous"){

    removeUI(selector = "div#milineagec_categorical_selector > div")
    removeUI(selector = "div#milineagec_subcategorical_selector > div")        

    insertUI(
      selector = "#milineagec_continuous_selector",
      where = "afterBegin",
      ui =  selectInput("milineagec_continuous", "Select a continuous", choices = sample_col_choices, 
        selected = NULL, multiple = FALSE)
    )
    
  }  else{

  }
 
})

observeEvent({physeq()}, { 
# this is a follow up of observeEvent({input$milineagec_arrange}
# when observeEvent({input$milineagec_arrange} first runs, it does not wait for physeq()to be not 
# NULL, so the output UI's are all empty. Thus, we need this observe event to wait untill physeq()
# finishes, and update the empty UI
 
  if(input$milineagec_arrange == "categorical"){ 

    updateSelectInput(session, "milineagec_categorical",
           choices = c("NULL", sample_variables(physeq())),
           selected = "NULL"
         )         
  }       
})

observeEvent({input$milineagec_id}, { 
# this observes milineagec_id and update timebase accordingly
  source("functions/mistudio_milineagec_timebase.R")
  updateSelectInput(session, "milineagec_timebase",
         choices = append(mistudio_milineagec_timebase(physeq(), input$milineagec_id), list("NULL" = "NULL"), 0),
         selected = "NULL"
  )        
})

observeEvent({input$milineagec_categorical}, { # observe what function user selects
      if(is.null(input$milineagec_categorical) || input$milineagec_categorical == "" || input$milineagec_categorical == "NULL"){
        updateSelectInput(session, "milineagec_subcategorical",
          choices = "",
          selected = NULL
        )  
        return()
      }

      paste_value_column = function(colname){
        sample = as.data.frame(sample_data(physeq())@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
          col.names = sample_data(physeq())@names, row.names = sample_data(physeq())@row.names) # manually reassign column and row names    
        return(paste0(sample[,colname], "%mistudio_seperator%", colname))
      }
      subcategorical_choices = unique(sapply(input$milineagec_categorical, paste_value_column)) 
      get_value = function(colname){
        sample = as.data.frame(sample_data(physeq())@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
          col.names = sample_data(physeq())@names, row.names = sample_data(physeq())@row.names) # manually reassign column and row names    
        return(sample[,colname])
      }

      names(subcategorical_choices) = unique(sapply(input$milineagec_categorical, get_value))
      updateSelectInput(session, "milineagec_subcategorical",
        choices = subcategorical_choices,
        selected = subcategorical_choices
      )    
     
})

observeEvent(ignoreNULL = FALSE,{input$milineagec_stratify}, { # observe what function user selects
      if(is.null(input$milineagec_stratify) || input$milineagec_stratify == ""){
        updateSelectInput(session, "milineagec_substratify",
          choices = "",
          selected = NULL
        )   
      }

      paste_value_column = function(colname){        
        sample = as.data.frame(sample_data(physeq())@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
          col.names = sample_data(physeq())@names, row.names = sample_data(physeq())@row.names) # manually reassign column and row names    
        return(paste0(sample[,colname], "%mistudio_seperator%", colname))
      }
      
      substratify_choices = unique(unlist(lapply(input$milineagec_stratify, paste_value_column)))
 
      if(!is.null(substratify_choices)){
        names(substratify_choices) = unlist(lapply(strsplit(substratify_choices, split = "%mistudio_seperator%"), '[', 1)) # get the first token of substratify_choices
      }
      

      updateSelectInput(session, "milineagec_substratify",
        choices = substratify_choices,
        selected = substratify_choices
      )    
     
})



# reactive .R functions



mistudio_milineagec_cova_result = reactive({

  source("functions/mistudio_milineagec_cova.R")

  mistudio_milineagec_cova(physeq())

})

mistudio_milineagec_conf_result = reactive({

  source("functions/mistudio_milineagec_conf.R")

  mistudio_milineagec_conf(physeq())

})

mistudio_milineagec_id_result = reactive({
  req(physeq())
  source("functions/mistudio_milineagec_id.R")

  mistudio_milineagec_id(physeq())

})
mistudio_milineagec_result = reactive({
  source("functions/mistudio_milineagec.R")

  if(input$action_milineagec == 0){ # return NULL if run button not clicked
    return(NULL)
  }

  isolate({

    if(input$milineagec_timebase == "NULL"){
      milineagec_id_parse = NULL
    } else{
      milineagec_id_parse = input$milineagec_timebase
    }

    if(input$milineagec_permtype == "NULL"){
      milineagec_permtype_parse = NULL
    } else {
      milineagec_permtype_parse = input$milineagec_permtype
    }
    result = mistudio_milineagec(
                phyobj = physeq(),
                milineagec_function = input$milineagec_function,
                milineagec_cova = input$milineagec_cova,
                milineagec_conf = input$milineagec_conf,
                milineagec_cova_extra = input$milineagec_cova_extra,
                milineagec_conf_extra = input$milineagec_conf_extra,
                milineagec_mindepth = input$milineagec_mindepth,
                milineagec_nresample = input$milineagec_nresample,
                milineagec_fdralpha = input$milineagec_fdralpha,
                milineagec_id = input$milineagec_id,
                milineagec_timebase = milineagec_id_parse, # could be "NULL", so extra parsing need to convert to NULL
                milineagec_permtype = milineagec_permtype_parse, # could be "NULL", so extra parsing need to convert to NULL
                milineagec_testtype  = input$milineagec_testtype
             )
    return(result)
  })
})
# mistudio_milineagec_compositional_result = eventReactive(input$action_milineagec_compositional,{
  
# })

mistudio_milineagec_compositional_result = eventReactive(input$action_milineagec_compositional,{
  source("functions/mistudio_milineagec_compositional.R")
  if(is.null(input$milineagec_categorical) || input$milineagec_categorical == "NULL"){ # it is possible at the beginning input$milienagec_continuous is "NULL" when physeq() not loaded, but it will become NULL once loaded. So we handle two situatinos here.
    categorical = NULL
  } else {
    categorical = input$milineagec_categorical
  }
  if(is.null(input$milineagec_continuous) || input$milineagec_continuous == "NULL"){ # it is possible at the beginning input$milienagec_continuous is "NULL" when physeq() not loaded, but it will become NULL once loaded. So we handle two situatinos here.
    continuous = NULL
  } else {
    continuous = input$milineagec_continuous
  }
  if(length(input$milineagec_stratify) > 2){
  	warning("more than two stratification found, using the first two")
  	stratify = input$milineagec_stratify[1:2]  	
  } else {
  	stratify = input$milineagec_stratify
  }
  if(input$milineagec_title == ""){
    title = NULL
  } else {
    title = input$milineagec_title
  }
  if(input$milineagec_xlab == ""){
    xlab = NULL
  } else {
    xlab = input$milineagec_xlab
  }
  if(input$milineagec_ylab == ""){
    ylab = NULL
  } else {
    ylab = input$milineagec_ylab
  }
  mistudio_milineagec_compositional(phyobj = physeq(), arrange = input$milineagec_arrange, lineage= input$milineagec_compositional_lineage, 
                                  stratify = stratify, substratify = input$milineagec_substratify, continuous = continuous,
                                  categorical = categorical, subcategorical = input$milineagec_subcategorical,
                                  proportion = input$milineagec_proportion,
                                  title = input$milineagec_title, xlab = xlab, ylab = ylab, palette = input$milineagec_palette, theme = input$milineagec_theme, hidex = input$milineagec_hidex
  )  
})


# UI outputs
output$milineagec_out_cova <- renderUI({
  req(physeq())
  selectInput("milineagec_cova", "Covariate(must select at least one covariate)", mistudio_milineagec_cova_result(), multiple = TRUE)
})

output$milineagec_out_conf <- renderUI({
  req(physeq())
  selectInput("milineagec_conf", "Confouding", mistudio_milineagec_conf_result(), multiple = TRUE)
})

output$milineagec_out_id <- renderUI({
  req(physeq())
  selectInput("milineagec_id", "ID", mistudio_milineagec_id_result(), selected = NULL, multiple = FALSE)
})

output$milineagec_out_timebase <- renderUI({
  req(physeq())
  selectInput("milineagec_timebase", "Time base\n do not select the whole category, select only sub categorical", c("NULL", mistudio_milineagec_cova_result()), selected = NULL, multiple = FALSE)
})

output$milineagec_out_permtype <- renderUI({
  selectInput("milineagec_permtype", "Perm Type", choices = c("NULL", "WTH", "BTW"), selected = NULL, multiple = FALSE)
})

output$milineagec_out_testtype <- renderUI({
  selectInput("milineagec_testtype", "Test Type", choices = c("chisq", "mix"), selected = NULL, multiple = FALSE)
})

output$milineagec_out_compositional_lineage <- renderUI({
	req(physeq())
  selectInput("milineagec_compositional_lineage", "Lineage(compositional plot)",choices = unique(na.omit(unlist(tax_table(physeq())))), multiple = FALSE)
})
output$milineagec_out_stratify <- renderUI({
  req(physeq())
  selectInput("milineagec_stratify", "Stratify", choices = sample_variables(physeq()), selected = NULL, multiple = TRUE)
})

output$milineagec_out_substratify <- renderUI({
  selectInput("milineagec_substratify", "Selected subcategories", choices = NULL, selected = NULL, multiple = TRUE)
})

# outputs
output$milineagec_lineagepval = renderTable(include.rownames=TRUE,
{
  get_df = function(QCAT_GEE_result){
    if(input$milineagec_function == "QCAT.Cluster"){

      return(as.data.frame(QCAT_GEE_result$lineage.pval, stringsAsFactors = FALSE))
    }
    rep_times = nrow(QCAT_GEE_result$lineage.pval$'Two-Part')
    parts = c("Two-Part", "Two-Part", "Zero-Part", "Zero-Part", "Postive-Part", "Postive-Part")
    types = c("Asymtopic", "Resampling","Asymtopic", "Resampling","Asymtopic", "Resampling")
    df_info = data.frame(Parts = parts, Types = types)

    mat = rbind(QCAT_GEE_result$lineage.pval$`Two-Part`, QCAT_GEE_result$lineage.pval$`Zero-Part`, QCAT_GEE_result$lineage.pval$`Positive-Part`)
    df_data = as.data.frame(mat, row.names = c('1','2','3','4','5','6'))

    result = try(cbind(df_info, df_data), silent = TRUE) # this line creates warnings when df_info or df_data
                                                       # are null, but it's expected. So warnings suppressed.
   
    return(result)
  }
  if(is.null(mistudio_milineagec_result())){
    return(NULL)
  }
  get_df(mistudio_milineagec_result())
}, digits = 10
)

output$milineagec_globalpval = renderTable(
{
  pvalue = mistudio_milineagec_result()$global.pval
  result = t(as.data.frame(pvalue, col.names = colnames(result$global.pval)))

}, digits = 10
)

output$milineagec_siglineage = renderTable({

  get_df = function(QCAT_GEE_result){
    if(is.null(QCAT_GEE_result$sig.lineage) || length(QCAT_GEE_result$sig.lineage) == 0){
      return("No significant lineage")
    }
    if(input$milineagec_function == "QCAT.Cluster"){

      return('Significant lineage' = as.data.frame(QCAT_GEE_result$sig.lineage))
    }
    length1 = length(QCAT_GEE_result$sig.lineage$`Two-Part`)
    length2 = length(QCAT_GEE_result$sig.lineage$`Zero-Part`)
    length3 = length(QCAT_GEE_result$sig.lineage$`Positive-Part`)

    model = character(max(length1, length2, length3))

    model1 = model
    model2 = model
    model3 = model

    if(length1 > 0) model1[1:length1] =  QCAT_GEE_result$sig.lineage$`Two-Part`
    if(length2 > 0) model2[1:length2] =  QCAT_GEE_result$sig.lineage$`Zero-Part`
    if(length3 > 0) model3[1:length3] =  QCAT_GEE_result$sig.lineage$`Positive-Part`

    result = data.frame('Two-Part' = model1, 'Zero-Part' = model2, 'Positive-Part' = model3)
    return(result)
  }
  get_df(mistudio_milineagec_result())
})
output$milineagec_compositional = renderPlot({
  mistudio_milineagec_compositional_result()
})
output$milineagec_compositional_reproduce <- downloadHandler(	
  filename = "mistudio_milineagec_reproduce",
  content = function(file) {
  	# save plot 
  	# ggsave(filename = "reproduce/reproduce_milineagec_compositional/milineagec_compositional.png", plot = mistudio_milineagec_compositional_result(), device = "png")
  	source("functions/mistudio_milineagec_compositional_reproduce.R")
  	# write files
  	if(is.null(input$milineagec_categorical) || input$milineagec_categorical == "NULL"){ # it is possible at the beginning input$milienagec_continuous is "NULL" when physeq() not loaded, but it will become NULL once loaded. So we handle two situatinos here.
  	  categorical = NULL
  	} else {
  	  categorical = input$milineagec_categorical
  	}
  	if(is.null(input$milineagec_continuous) || input$milineagec_continuous == "NULL"){ # it is possible at the beginning input$milienagec_continuous is "NULL" when physeq() not loaded, but it will become NULL once loaded. So we handle two situatinos here.
  	  continuous = NULL
  	} else {
  	  continuous = input$milineagec_continuous
  	}
  	if(input$milineagec_xlab == ""){
  	  xlab = NULL
  	} else {
  	  xlab = input$milineagec_xlab
  	}
  	if(input$milineagec_ylab == ""){
  	  ylab = NULL
  	} else {
  	  ylab = input$milineagec_ylab  
  	}
  	mistudio_milineagec_compositional_reproduce(
  			phyobj_original = get_phyloseq_data(),
  		    phyobj_filtered = physeq(),
  		    arrange = input$milineagec_arrange, 
  		    lineage= input$milineagec_compositional_lineage, 
  		    stratify = input$milineagec_stratify, 
  		    substratify = input$milineagec_substratify, 
  		    continuous = continuous,
  		    categorical = categorical, 
  		    subcategorical = input$milineagec_subcategorical,
  		    xlab = xlab, 
  		    ylab = ylab, 
  		    palette = input$milineagec_palette, 
  		    theme = input$milineagec_theme, 
  		    hidex = input$milineagec_hidex
  	)
    file.copy("reproduce/reproduce_milineagec_compositional/mistudio_milineagec_compositional_reproduce.zip", file) # in this and the next line, we need to move the file to the path indicated by "file" in order for downloadHandler to work
    file.remove("reproduce/reproduce_milineagec_compositional/mistudio_milineagec_compositional_reproduce.zip")
  },
  contentType = "application/zip"
)
outputOptions(output, "milineagec_compositional", suspendWhenHidden = FALSE)  
