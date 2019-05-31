# make dynamically slection menu based on functions being selected
observeEvent({input$milineage_function}, { # observe what function user selects
      
      if(input$milineage_function == "QCAT_GEE"){ 
        updateSelectInput(session, "milineage_cova",
              label = "Covariate for positive part"
              )

        updateSelectInput(session, "milineage_conf",
              label = "Confouding for positive part"
              )
        
        insertUI(
          selector = "#milineage_cova_selector",
          where = "afterBegin",
          ui =  selectInput("milineage_cova_extra", "Covariate for zero part", mistudio_milineage_cova_result(), 
            selected = input$milineage_cova, multiple = TRUE)
        )  

        insertUI(
          selector = "#milineage_conf_selector",
          where = "afterBegin",
          ui =  selectInput("milineage_conf_extra", "Confounding for zero part", mistudio_milineage_conf_result(), 
            selected = input$milineage_conf, multiple = TRUE)
        )         
        removeUI(selector = "div#milineage_ZILB_selector > div")
        removeUI(selector = "div#milineage_testtype_selector > div") 
        
      } else if (input$milineage_function == "ZIGDM"){
        updateSelectInput(session, "milineage_cova",
              label = "Covariate"
              )

        updateSelectInput(session, "milineage_conf",
              label = "Confouding"
              )
        removeUI(selector = "div#milineage_cova_selector > div")
        removeUI(selector = "div#milineage_conf_selector > div")        

        insertUI(
          selector = "#milineage_ZILB_selector",
          where = "afterBegin",
          ui =  numericInput("milineage_ZI.LB", "Zero-inflated model lower bound", 10)          
        )

        insertUI(
          selector = "#milineage_testtype_selector",
          where = "afterBegin",
          ui =  selectInput("milineage_testtype", "Test Type", choices = c("mean abundance", "dispersion level", "presence-absence frequency"), multiple = FALSE)
        ) 
        
      } else if(input$milineage_function == "QCAT"){
        updateSelectInput(session, "milineage_cova",
              label = "Covariate"
              )

        updateSelectInput(session, "milineage_conf",
              label = "Confouding"
              )
        removeUI(selector = "div#milineage_cova_selector > div")
        removeUI(selector = "div#milineage_conf_selector > div")
        removeUI(selector = "div#milineage_ZILB_selector > div")
        removeUI(selector = "div#milineage_testtype_selector > div")

      } else{

      }



      
    })

# let cova(conf) and cova_extra(conf_extra) always have the same value unless explicily modified
observeEvent(c(input$milineage_cova, input$milineage_conf), { # using c() is the correct way to observe two events
  # req(input$milineage_cova_extra)
  # req(input$milineage_conf_extra)
  updateSelectInput(session, "milineage_cova_extra",
        selected = input$milineage_cova
        )

  updateSelectInput(session, "milineage_conf_extra",
        selected = input$milineage_conf
        )
})

observeEvent({input$milineage_arrange}, { # observe what function user selects
  if(is.null(physeq())){
    sample_col_choices = NULL
  } else {
    sample_col_choices = c("NULL", sample_variables(physeq()))
  }

  if(input$milineage_arrange == "categorical"){ 

    removeUI(selector = "div#milineage_continuous_selector > div")
    
  
    insertUI(
      selector = "#milineage_categorical_selector",
      where = "afterBegin",
      ui =  selectInput("milineage_categorical", "Group and average sample by", choices = sample_col_choices, 
        selected = NULL, multiple = FALSE)
    ) 
    insertUI(
      selector = "#milineage_subcategorical_selector",
      where = "afterBegin",
      ui =  selectInput("milineage_subcategorical", "Selected categories", choices = NULL, 
        selected = NULL, multiple = TRUE)
    )  

    
  } else if (input$milineage_arrange == "continuous"){

    removeUI(selector = "div#milineage_categorical_selector > div")
    removeUI(selector = "div#milineage_subcategorical_selector > div")        

    insertUI(
      selector = "#milineage_continuous_selector",
      where = "afterBegin",
      ui =  selectInput("milineage_continuous", "Select a continuous", choices = sample_col_choices, 
        selected = NULL, multiple = FALSE)
    )
    
  }  else{

  }
 
})

observeEvent({physeq()}, { 
# this is a follow up of observeEvent({input$milineage_arrange}
# when observeEvent({input$milineage_arrange} first runs, it does not wait for physeq()to be not 
# NULL, so the output UI's are all empty. Thus, we need this observe event to wait untill physeq()
# finishes, and update the empty UI
 
  if(input$milineage_arrange == "categorical"){ 

    updateSelectInput(session, "milineage_categorical",
           choices = c("NULL", sample_variables(physeq())),
           selected = "NULL"
         )         
  }       
})

observeEvent({input$milineage_categorical}, { # observe what function user selects
      if(is.null(input$milineage_categorical) || input$milineage_categorical == "" || input$milineage_categorical == "NULL"){
        updateSelectInput(session, "milineage_subcategorical",
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
      subcategorical_choices = unique(sapply(input$milineage_categorical, paste_value_column)) 
      get_value = function(colname){
        sample = as.data.frame(sample_data(physeq())@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
          col.names = sample_data(physeq())@names, row.names = sample_data(physeq())@row.names) # manually reassign column and row names    
        return(sample[,colname])
      }

      names(subcategorical_choices) = unique(sapply(input$milineage_categorical, get_value))
      updateSelectInput(session, "milineage_subcategorical",
        choices = subcategorical_choices,
        selected = subcategorical_choices
      )    
     
})

observeEvent(ignoreNULL = FALSE,{input$milineage_stratify}, { # observe what function user selects
      if(is.null(input$milineage_stratify) || input$milineage_stratify == ""){
        updateSelectInput(session, "milineage_substratify",
          choices = "",
          selected = NULL
        )   
      }

      paste_value_column = function(colname){        
        sample = as.data.frame(sample_data(physeq())@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
          col.names = sample_data(physeq())@names, row.names = sample_data(physeq())@row.names) # manually reassign column and row names    
        return(paste0(sample[,colname], "%mistudio_seperator%", colname))
      }
      
      substratify_choices = unique(unlist(lapply(input$milineage_stratify, paste_value_column)))
 
      if(!is.null(substratify_choices)){
        names(substratify_choices) = unlist(lapply(strsplit(substratify_choices, split = "%mistudio_seperator%"), '[', 1)) # get the first token of substratify_choices
      }
      

      updateSelectInput(session, "milineage_substratify",
        choices = substratify_choices,
        selected = substratify_choices
      )    
     
})


# reactive .R functions
mistudio_milineage_cova_result = reactive({
  req(physeq())

  source("functions/mistudio_milineage_cova.R")

  mistudio_milineage_cova(physeq())

})

mistudio_milineage_conf_result = reactive({
  req(physeq())

  source("functions/mistudio_milineage_conf.R")

  mistudio_milineage_conf(physeq())

})

mistudio_milineage_result = reactive({

  source("functions/mistudio_milineage.R")

  if(input$action_milineage == 0){ # return NULL if run button not clicked
    return(NULL)
  }
  isolate({
    result = mistudio_milineage(
        phyobj = physeq(),
        milineage_function = input$milineage_function,
        milineage_cova = input$milineage_cova,
        milineage_conf = input$milineage_conf,
        milineage_cova_extra = input$milineage_cova_extra,
        milineage_conf_extra = input$milineage_conf_extra,
        milineage_mindepth = input$milineage_mindepth,
        milineage_nresample = input$milineage_nresample,
        milineage_fdralpha = input$milineage_fdralpha,
        milineage_ZI.LB = input$milineage_ZI.LB,
        milineage_testtype  = input$milineage_testtype
    )

    return(result)
  })
})

mistudio_milineage_circular_result = eventReactive(input$action_milineage_compositional,{
  source("functions/mistudio_milineage_circular.R")
  mistudio_milineage_circular(phyobj = physeq(), lineage = input$milineage_annot_lineage)  
})

mistudio_milineage_compositional_result = eventReactive(input$action_milineage_compositional,{
  source("functions/mistudio_milineage_compositional.R")
  if(is.null(input$milineage_categorical) || input$milineage_categorical == "NULL"){ # it is possible at the beginning input$milienagec_continuous is "NULL" when physeq() not loaded, but it will become NULL once loaded. So we handle two situatinos here.
    categorical = NULL
  } else {
    categorical = input$milineage_categorical
  }
  if(is.null(input$milineage_continuous) || input$milineage_continuous == "NULL"){ # it is possible at the beginning input$milienagec_continuous is "NULL" when physeq() not loaded, but it will become NULL once loaded. So we handle two situatinos here.
    continuous = NULL
  } else {
    continuous = input$milineage_continuous
  }
  if(length(input$milineage_stratify) > 2){
  	warning("more than two stratification found, using the first two")
  	stratify = input$milineage_stratify[1:2]  	
  } else {
  	stratify = input$milineage_stratify
  }
  if(input$milineage_title == ""){
    title = NULL
  } else {
    title = input$milineage_title
  }
  if(input$milineage_xlab == ""){
    xlab = NULL
  } else {
    xlab = input$milineage_xlab
  }
  if(input$milineage_ylab == ""){
    ylab = NULL
  } else {
    ylab = input$milineage_ylab
  }

  mistudio_milineage_compositional(phyobj = physeq(), arrange = input$milineage_arrange, lineage= input$milineage_compositional_lineage, 
                                  stratify = stratify, substratify = input$milineage_substratify, continuous = continuous,
                                  categorical = categorical, subcategorical = input$milineage_subcategorical,
                                  proportion = input$milineage_proportion,
                                  title = input$milineage_title, xlab = xlab, ylab = ylab, palette = input$milineage_palette, theme = input$milineage_theme, hidex = input$milineage_hidex
  )  
})

# UI outputs
output$milineage_out_cova <- renderUI({
  selectInput("milineage_cova", "Covariate", mistudio_milineage_cova_result(), multiple = TRUE)
})

output$milineage_out_conf <- renderUI({
  selectInput("milineage_conf", "Confouding", mistudio_milineage_conf_result(), multiple = TRUE)
})
output$milineage_out_annot_lineage <- renderUI({
  selectInput("milineage_annot_lineage", "Lineage(tree plot)", choices = unique(na.omit(unlist(tax_table(physeq())))), multiple = TRUE)

})
output$milineage_out_compositional_lineage <- renderUI({
  selectInput("milineage_compositional_lineage", "Lineage(compositional plot)",choices = unique(na.omit(unlist(tax_table(physeq())))), multiple = TRUE)
})
output$milineage_out_stratify <- renderUI({
  req(physeq())
  selectInput("milineage_stratify", "Stratify", choices = sample_variables(physeq()), selected = NULL, multiple = TRUE)
})

output$milineage_out_substratify <- renderUI({
  selectInput("milineage_substratify", "Selected subcategories", choices = NULL, selected = NULL, multiple = TRUE)
})




# outputs
output$lineage.pval = renderTable(
{
  get_df = function(QCAT_GEE_result){
    if(input$milineage_function == "QCAT" || input$milineage_function == "ZIGDM"){
      return(as.data.frame(QCAT_GEE_result$lineage.pval, stringsAsFactors = FALSE))
    }
    # in case of user first perform QCAT, then select QCAT_GEE, the p-value display will have error 
    # because 'Two-Part' does not exist in the result. Thus, we return NULL when it happens
    if(class(try(QCAT_GEE_result$lineage.pval$'Two-Part', silent = TRUE)) == "try-error"){
      return(NULL)
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
  if(is.null(mistudio_milineage_result())){
  	return(NULL)
  }
  get_df(mistudio_milineage_result())
}, digits = 10
  # matrix(c(as.character(names(mistudio_milineage_result()$global.pval)), as.character(mistudio_milineage_result()$global.pval)), byrow = T, nrow = 2)
)

output$global.pval = renderTable(
{
	if(is.null(mistudio_milineage_result())){
		return(NULL)
	}
  pvalue = mistudio_milineage_result()$global.pval
  result = t(as.data.frame(pvalue, col.names = colnames(result$global.pval)))

}, digits = 10
  # matrix(c(as.character(names(mistudio_milineage_result()$global.pval)), as.character(mistudio_milineage_result()$global.pval)), byrow = T, nrow = 2)
)

output$sig.lineage = renderTable({

  get_df = function(QCAT_GEE_result){
    if(is.null(QCAT_GEE_result$sig.lineage) || length(QCAT_GEE_result$sig.lineage) == 0){
      return("No significant lineage")
    }

    if(input$milineage_function == "QCAT" || input$milineage_function == "ZIGDM"){

      updateSelectInput(session, "milineage_compositional_lineage",
                    selected = QCAT_GEE_result$sig.lineage[1]
                    )
      updateSelectInput(session, "milineage_annot_lineage",
                    selected = QCAT_GEE_result$sig.lineage
                    )
      
      return('Significant lineage' = as.data.frame(QCAT_GEE_result$sig.lineage))
    }
    updateSelectInput(session, "milineage_compositional_lineage",
                  selected = QCAT_GEE_result$sig.lineage[1]
                  )
    updateSelectInput(session, "milineage_annot_lineage",
                  selected = QCAT_GEE_result$sig.lineage
                  )    
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
  if(is.null(mistudio_milineage_result())){
  	return(NULL)
  }
  get_df(mistudio_milineage_result())
  # matrix(c(as.character(names(mistudio_milineage_result()$sig.lineage)), as.character(mistudio_milineage_result()$sig.lineage)), byrow = T, nrow = 2)
})


output$milineage_tree <- renderPlot({
  mistudio_milineage_circular_result()    
})


output$milineage_compositional = renderPlot({
  result = mistudio_milineage_compositional_result()
  return(result)
})


# Download Provenence 

output$milineage_reproduce <- downloadHandler(	
  filename = "mistudio_milineage_reproduce",
  content = function(file) {
  	source("functions/mistudio_milineage_reproduce.R")

  	# clean up input
    filter_rank_selection = av(input$filter_rank_selection)
    filter_rank = av(input$filter_rank)
    filter_samvars_selection = av(input$filter_samvars_selection)
    filter_samvars = av(input$filter_samvars)
    if(is.null(input$milineage_categorical) || input$milineage_categorical == "NULL"){ # it is possible at the beginning input$milienagec_continuous is "NULL" when physeq() not loaded, but it will become NULL once loaded. So we handle two situatinos here.
      categorical = NULL
    } else {
      categorical = input$milineage_categorical
    }
    if(is.null(input$milineage_continuous) || input$milineage_continuous == "NULL"){ # it is possible at the beginning input$milienagec_continuous is "NULL" when physeq() not loaded, but it will become NULL once loaded. So we handle two situatinos here.
      continuous = NULL
    } else {
      continuous = input$milineage_continuous
    }
    if(length(input$milineage_stratify) > 2){
      warning("more than two stratification found, using the first two")
      stratify = input$milineage_stratify[1:2]   
    } else {
      stratify = input$milineage_stratify
    }
    if(input$milineage_xlab == ""){
      xlab = NULL
    } else {
      xlab = input$milineage_xlab
    }
    if(input$milineage_ylab == ""){
      ylab = NULL
    } else {
      ylab = input$milineage_ylab
    }

    mistudio_milineage_reproduce(
      path = "reproduce/milineage/milineage_reproduce.R",
      phyobj_name = input$physeqSelect,
      data_path = "data",

      filter_rank_selection,
      filter_rank,
      filter_samvars_selection,
      filter_samvars,
      input$filter_otu_aggregate,
      input$filter_otu_zerofraction,
      input$filter_otu_minsum,

      input$milineage_function,
      input$milineage_cova,
      input$milineage_conf,
      input$milineage_cova_extra,
      input$milineage_conf_extra,
      input$milineage_mindepth,
      input$milineage_nresample,
      input$milineage_fdralpha,
      input$milineage_ZI.LB,
      input$milineage_testtype,

      arrange = input$milineage_arrange, 
      lineage= input$milineage_compositional_lineage,     
      stratify = stratify, 
      substratify = input$milineage_substratify, 
      continuous = continuous,      
      categorical = categorical, 
      subcategorical = input$milineage_subcategorical,
      xlab = xlab, 
      ylab = ylab, 
      palette = input$milineage_palette, 
      theme = input$milineage_theme, 
      hidex = input$milineage_hidex
      )
    	# mistudio_milineage_reproduce(
    	# 		phyobj_original = get_phyloseq_data(),
    	# 	    phyobj_filtered = physeq(),
    	# 	    milineage_function = input$milineage_function,
    	# 	    milineage_cova = input$milineage_cova,
    	# 	    milineage_conf = input$milineage_conf,
    	# 	    #milineage_cova_extra = milineage_cova_extra,
    	# 	    milineage_cova_extra = input$milineage_cova_extra,
    	# 	    #milineage_conf_extra = milineage_conf_extra,
    	# 	    milineage_conf_extra = input$milineage_conf_extra,
    	# 	    milineage_mindepth = input$milineage_mindepth,
    	# 	    milineage_nresample = input$milineage_nresample,
    	# 	    milineage_fdralpha = input$milineage_fdralpha,
    	# 	    milineage_ZI.LB = input$milineage_ZI.LB,
    	# 	    milineage_testtype = input$milineage_testtype)
    file.copy(paste0(input$physeqSelect, "_milineage", "_reproduce", ".zip"), file) # in this and the next line, we need to move the file to the path indicated by "file" in order for downloadHandler to work
    file.remove(paste0(input$physeqSelect, "_milineage", "_reproduce", ".zip"))
  },
  contentType = "application/zip"
)
outputOptions(output, "milineage_tree", suspendWhenHidden = FALSE)  
outputOptions(output, "milineage_compositional", suspendWhenHidden = FALSE)  
