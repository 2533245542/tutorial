# reactive UI
# library(shiny)

COVA_COL = reactive({
  as.list(sample_data(physeq())@names)
})

COVA = reactive({
  as.list(sample_data(physeq()))
})

mistudio_milineage_result = reactive({

  source("functions/mistudio_milineage.R")

  if(input$action_milineage == 0){ # return NULL if run button unclicked
    return(NULL)
  }

  isolate({
    result = mistudio_milineage(physeq(), input$milineage_cov_n, input$milineage_cov_c, input$milineage_con_n, 
      input$milineage_con_c, input$milineage_nresample, input$milineage_fdralpha, input$milineage_mindepth)

    return(result)
  })
})

mistudio_milineage_compositional_result = reactive({
  source("functions/mistudio_milineage_compositional.R")

  if(input$action_milineage_compositional == 0){
    return(NULL)
  }

  isolate({

    # stopifnot(!is.null(input$milineage_lineage))
    # stopifnot(!is.null(mistudio_milineage_result()))

    result = mistudio_milineage_compositional(phyobj = physeq(), lineage = input$milineage_lineage, sortby = input$milineage_sortby, 
      stratifyby = input$milineage_stratifyby)
    # print(str(result))
    return(result)
  })

})





output$milineage_out_cov_n <- renderUI({
  selectInput("milineage_cov_n", "Numeric Covariate", COVA_COL(), multiple = TRUE)
})

output$milineage_out_cov_c <- renderUI({
  selectInput("milineage_cov_c", "Categorical Covariate", COVA(), multiple = TRUE)
})

output$milineage_out_con_n <- renderUI({
  selectInput("milineage_con_n", "Numeric Confouding", COVA_COL(), multiple = TRUE)
})

output$milineage_out_con_c <- renderUI({
  selectInput("milineage_con_c", "Categorical Confouding", COVA(), multiple = TRUE)
})

output$milineage_out_lineage <- renderUI({
  selectInput("milineage_lineage", "Select Lineage To Plot", choices = unique(as.character(tax_table(physeq()))))
})

output$milineage_out_sortby <- renderUI({
  choices = append(COVA_COL(), c("NULL" = "NULL"), 0) # vector name must have at least a white space 
  return(selectInput("milineage_sortby", "X Axis Sorted By", choices = choices))
})

output$milineage_out_stratifyby <- renderUI({
  choices = append(COVA_COL(), c("NULL" = "NULL"), 0)
  return(selectInput("milineage_stratifyby", "Sample Stratify By", choices = choices))
})






output$lineage.pval = renderTable(
{
  get_df = function(QCAT_GEE_result){
    rep_times = nrow(QCAT_GEE_result$lineage.pval$'Two-Part')
    parts = c("Two-Part", "Two-Part", "Zero-Part", "Zero-Part", "Postive-Part", "Postive-Part")
    types = c("Asymtopic", "Resampling","Asymtopic", "Resampling","Asymtopic", "Resampling")
    df_info = data.frame(Parts = parts, Types = types)

    mat = rbind(QCAT_GEE_result$lineage.pval$`Two-Part`, QCAT_GEE_result$lineage.pval$`Zero-Part`, QCAT_GEE_result$lineage.pval$`Positive-Part`)
    df_data = as.data.frame(mat, row.names = c('1','2','3','4','5','6'))

    suppressWarnings(result = cbind(df_info, df_data)) # this line creates warnings when df_info or df_data
                                                       # are null, but it's expected. So warnings suppressed.
    return(result)
  }
  get_df(mistudio_milineage_result())
}, digits = 10
  # matrix(c(as.character(names(mistudio_milineage_result()$global.pval)), as.character(mistudio_milineage_result()$global.pval)), byrow = T, nrow = 2)
)

output$global.pval = renderTable(
{
  pvalue = mistudio_milineage_result()$global.pval
  result = t(as.data.frame(pvalue, col.names = colnames(result$global.pval)))

}, digits = 10
  # matrix(c(as.character(names(mistudio_milineage_result()$global.pval)), as.character(mistudio_milineage_result()$global.pval)), byrow = T, nrow = 2)
)

output$sig.lineage = renderTable({
  get_df = function(QCAT_GEE_result){
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
  get_df(mistudio_milineage_result())
  # matrix(c(as.character(names(mistudio_milineage_result()$sig.lineage)), as.character(mistudio_milineage_result()$sig.lineage)), byrow = T, nrow = 2)
})

output$myImage <- renderImage({
  mistudio_milineage_result()
  filename = ""
  list(src = filename)
}, deleteFile = FALSE)


output$milineage_compositional = renderPlot({
  result = mistudio_milineage_compositional_result()
  return(result)
})


# Download Provenence 

fname = "reproduce_milineage.zip"
output$downloadData <- downloadHandler(
  filename <- function() {
    fname
  },
  content <- function(file) {
    # Generate records
    source("functions_milineage_repro_script.R")
    phyloseqobject = physeq()
    provenence_milineage(
                         input$milineage_cova, input$milineage_cova_col,
                         input$milineage_conf, input$milineage_conf_col,
                         physeq()
                         )
    # Generate zip file 
    fs = dir(path = "provenence_milineage/")
    fs = paste0("provenence_milineage/", fs)
    fname = paste0("provenence_milineage/", fname)
    zip(fname, fs)
    
    # Output zip file 
    file.copy(fname, file)
    
    # Remove zip file
    file.remove(fname)
  },
  contentType = "application/zip"
)
