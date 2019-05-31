# observeEvent({input$filter_subset_rank}, { 
#   if(!is.null(input$filter_subset_rank)){
#     taxa_list = lapply(input$filter_subset_rank, function(x) return(as.character(tax_table(get_phyloseq_data()) [, x])))
#     updateSelectInput(session, "filter_subset_taxa",
#           choices = unique(do.call(c, taxa_list))
#           )
#   }
# })
################################################################################
# UI subset_taxa expression cascade
# filter_subset_taxa_expr
################################################################################
output$filter_out_otu_aggregate <- renderUI({
  rankNames = list("NULL"="NULL")
  rankNames <- c(rankNames, as.list(rank_names(get_phyloseq_data(), errorIfNULL=FALSE)))
  rankNames <- c(rankNames, list(OTU="OTU"))
  return(
    selectInput("filter_otu_aggregate", "Aggregate", rankNames, "NULL", multiple = FALSE)
  )
})

output$filter_out_subset_rank <- renderUI({
  return(
    selectInput("filter_subset_rank", "Subset taxonomy table(based on rank)", choices = rank_names(get_phyloseq_data()), multiple = TRUE)
  )
})
output$filter_out_subset_taxa <- renderUI({  
  choices = unique(as.character(tax_table(get_phyloseq_data())))
  return(
    selectInput("filter_subset_taxa", "Subset taxonomy table(based on lineage)", choices =  choices, multiple = TRUE)
  )
})
output$filter_out_subset_subject <- renderUI({ # this one allows the user to select rows in the sample data to preserve
  return(
    selectInput("filter_subset_subject", "Subset meta data table(based on subjects)", choices = sample_names(get_phyloseq_data()), multiple = TRUE)
  )
})
output$filter_out_subset_variable <- renderUI({ # this one allows the user to select categorical variables presented in the sample table
  get_sample_list = function(phyobj){
    sample = as.data.frame(sample_data(phyobj)@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
      col.names = sample_data(phyobj)@names, row.names = sample_data(phyobj)@row.names) # manually reassign column and row names    
    names = colnames(sample)
    myfun = function(list, name){
      list_unique = unique(list)
      if(length(list_unique) != nrow(sample)){ # only make a list when not all values in a column is unique(is a categorical column)
        content = paste0(list_unique, "%mistudio_seperator%", name)
        names(content) = list_unique
        return(content)
      } else {
        return(NULL)
      }  
    }
    mapply(myfun, sample, names, SIMPLIFY = FALSE)
  }
  return(
    selectInput("filter_subset_variable", "Subset meta data table(based on variables)", choices =  get_sample_list(get_phyloseq_data()), multiple = TRUE)
  )
})

output$filter_uix_subset_taxa_ranks <- renderUI({
  rankNames = list("NULL"="NULL")
  rankNames <- c(rankNames, as.list(rank_names(get_phyloseq_data(), errorIfNULL=FALSE)))
  # rankNames <- c(rankNames, list(OTU="OTU"))
  return(
    selectInput("filter_rank", "Subset Ranks", rankNames, "NULL", multiple = FALSE)
  )
})
output$filter_uix_subset_taxa_select <- renderUI({
  rank_list = list("NULL" = "NULL")
  if(!is.null(av(input$filter_rank))){
    # If a filter rank is specified, use it, and provide the multi-select widget for these rank classes
    if(input$filter_rank == "OTU"){
      rank_list <- c(rank_list, as.list(taxa_names(get_phyloseq_data())))
    } else {
      rank_list <- c(rank_list, as.list(get_taxa_unique(get_phyloseq_data(), input$filter_rank)))
    }
  }
  return(
    selectInput(inputId = "filter_rank_selection", label = "Subset Taxa",
                choices = rank_list, selected = "NULL", multiple = TRUE)
  )
})
################################################################################
# UI subset_samples expression cascade
# filter_subset_samp_expr
################################################################################
output$filter_uix_subset_sample_vars <- renderUI({
  sampVars = list("NULL"="NULL")
  sampVars <- c(sampVars, as.list(sample_variables(get_phyloseq_data(), errorIfNULL=FALSE)))
  sampVars <- c(sampVars, list(Sample="Sample"))
  return(
    selectInput("filter_samvars", "Subset sample variables", sampVars, "NULL", multiple = FALSE)
  )
})
output$filter_uix_subset_sample_select <- renderUI({
  varLevels = list("NULL"="NULL")
  if(!is.null(av(input$filter_samvars))){
    if(input$filter_samvars == "Sample"){
      varLevels <- c(varLevels, as.list(sample_names(get_phyloseq_data())))
    } else {
      if(!is.null(sample_variables(get_phyloseq_data(), FALSE))){
        varvec = get_variable(get_phyloseq_data(), input$filter_samvars)
        if(plyr::is.discrete(varvec)){
          varLevels <- c(varLevels, as.list(unique(as(varvec, "character"))))
        }
      } 
    }
  }
  return(
    selectInput(inputId = "filter_samvars_selection", label = "Subset variable classes",
                choices = varLevels, selected = "NULL", multiple = TRUE)
  )  
})


output$filtered_contents0 <- renderUI({
  output_phyloseq_print_html(get_phyloseq_data())
})
output$filtered_contents <- renderUI({
  output_phyloseq_print_html(physeq())
})

# ################################################################################
# # Component Table
# ################################################################################
output$uix_available_components_filt <- renderUI({
  selectInput("available_components_filt", "Available Components",
              choices = component_options(physeq()))
})
# # Render the user-selected data component using DataTables
output$physeqComponentTable <- renderDataTable({
  if(is.null(av(input$available_components_filt))){
    return(NULL)
  }
  component = do.call(what = input$available_components_filt, args = list(physeq()))
  return(tablify_phyloseq_component(component, input$component_table_colmax_filt))
}, options = list(
  pageLength = 5 
))
################################################################################
# The main reactive data object. Returns a phyloseq-class instance.
# This is considered the "filtered" data, used by all downstream panels,
# And generally the input to any transformation options as well
################################################################################

physeq = reactive({
  req(get_phyloseq_data())
  ps0 = get_phyloseq_data()

  if(input$action_filter == 0){
    # Don't execute filter if filter-button has never been clicked.
    if(inherits(ps0, "phyloseq")){

      # the most important step is to remove factors in the dataframe of metadata
      drop_df_factor = function(df){
        char_cols_boolean = lapply(df, function(x) suppressWarnings(all(is.na(as.numeric(as.character(x))))))
        df_in_list = mapply(function(x, y){
          if(x){
            return(as.character(y))
          } else {
            return(as.numeric(y))
          }}, char_cols_boolean,df,SIMPLIFY = FALSE
        )
        return(as.data.frame(df_in_list, stringsAsFactors=FALSE, row.names = rownames(df)))
      }
      ps0 = phyloseq(sample_data(drop_df_factor(sample_data(ps0))), tax_table(ps0), otu_table(ps0))
      return(ps0)
    } else {
      return(NULL)
    }
  }
  # Isolate all filter code so that button click is required for update
  isolate({
    if(inherits(ps0, "phyloseq")){
      # the most important step is to remove factors in the dataframe of metadata
      drop_df_factor = function(df){
        char_cols_boolean = lapply(df, function(x) all(is.na(as.numeric(as.character(x)))))
        df_in_list = mapply(function(x, y){
          if(x){
            return(as.character(y))
          } else {
            return(as.numeric(y))
          }}, char_cols_boolean,df,SIMPLIFY = FALSE
        )
        return(as.data.frame(df_in_list, stringsAsFactors=FALSE, row.names = rownames(df)))
      }
      ps0 = phyloseq(sample_data(drop_df_factor(sample_data(ps0))), tax_table(ps0), otu_table(ps0))
      # Cascading selection filters
      if( !is.null(av(input$filter_rank_selection)) ){
        keepTaxa = NULL
        if(!is.null(tax_table(ps0, FALSE))){
          if(input$filter_rank == "OTU"){
            # OTU IDs directly
            keepTaxa = input$filter_rank_selection
          } else {
            TT = as(tax_table(ps0), "matrix")
            keepTaxa = TT[, input$filter_rank] %in% input$filter_rank_selection 
          }
          if(length(keepTaxa) > 1){
            ps0 <- prune_taxa(keepTaxa, ps0)
          } else {
            warning("Bad subset_taxa specification. ntaxa(ps0) one or fewer OTUs")
          }
        }
      }
      if(!is.null(av(input$filter_samvars_selection)) ){
        keepSamples = NULL
        if(!is.null(sample_data(ps0, FALSE))){
          if(input$filter_samvars == "Sample"){
            # Samples IDs directly
            keepSamples = input$filter_samvars_selection
          } else {
            varvec = as(get_variable(ps0, input$filter_samvars), "character")
            keepSamples = varvec %in% input$filter_samvars_selection 
          }
          if(length(keepSamples) > 1){
            ps0 <- prune_samples(keepSamples, ps0)
          } else {
            warning("Bad subset_taxa specification. ntaxa(ps0) one or fewer OTUs")
          }
        }
      }

      if((input$filter_otu_aggregate) != "NULL"){
        source("functions/mistudio_filter_otu_aggregate.R")
        ps0 = mistudio_filter_otu_aggregate(ps0, input$filter_otu_aggregate)
      }
      if(!is.null(input$filter_subset_rank) || !is.null(input$filter_subset_taxa)){
        source("functions/mistudio_filter_tax_subset.R")
        ps0 = mistudio_filter_tax_subset(ps0, input$filter_subset_rank, input$filter_subset_taxa)
      }
      if(!is.null(input$filter_subset_subject)){
        source("functions/mistudio_filter_sample_subset_subject.R")
        ps0 = mistudio_filter_sample_subset_subject(ps0, input$filter_subset_subject)
      }
      if(!is.null(input$filter_subset_variable)){
        source("functions/mistudio_filter_sample_subset_variable.R")
        ps0 = mistudio_filter_sample_subset_variable(ps0, input$filter_subset_variable)
      }
      if( input$filter_otu_minsum > 0 ){
        source("functions/mistudio_filter_otu_minsum.R")
        ps0 = mistudio_filter_otu_minsum(ps0, input$filter_otu_minsum)
      }
      if(input$filter_otu_zerofraction > 0 ){
        source("functions/mistudio_filter_otu_zerofraction.R")
        ps0 = mistudio_filter_otu_zerofraction(ps0, input$filter_otu_zerofraction)
      }

      return(ps0)
    } else {
      return(NULL)
    }
  })
})

output$action_export <- downloadHandler(   
  filename = "filtered_phyloseq_object.RData",
  content = function(file) {    
    req(physeq())
    filtered_phyloseq_obj = physeq()
    save(filtered_phyloseq_obj, file = "temp/filtered_phyloseq_obj.RData")    
    file.copy("temp/filtered_phyloseq_obj.RData", file) # in this and the next line, we need to move the file to the path indicated by "file" in order for downloadHandler to work
    file.remove("temp/filtered_phyloseq_obj.RData")
  },
  contentType = NULL
)

output$filter_plot_min_sum_before_filter = renderPlot({
  req(get_phyloseq_data())
  get_otu = function(phyobj){
    otu = otu_table(phyobj)@.Data
    tax = tax_table(phyobj)@.Data
    if(all(rownames(tax) == rownames(otu))){
      return(otu)
    } else {
      return(t(otu))
    }     
  }
  otu = t(get_otu(get_phyloseq_data()))
  return(
    qplot(rowSums(otu), geom = "histogram", 
      xlab = "Adding each sample's otu count together", 
      ylab = "The number of sample with that amount of count"
    )
  )
})

output$filter_plot_zerofraction_before_filter = renderPlot({
  req(get_phyloseq_data())
  get_otu = function(phyobj){
        otu = otu_table(phyobj)@.Data
        tax = tax_table(phyobj)@.Data     
        if(all(rownames(tax) == rownames(otu))){
          return(otu)
        } else {
          return(t(otu))
        }     
      }
  otu = t(get_otu(get_phyloseq_data()))
  otu_zerofraction = apply(otu, MARGIN = 2, function(x) return((sum(x == 0) / length(x))))       
  
  return(
    qplot(otu_zerofraction, geom = "histogram", 
      xlab = "Proportion of zero in a otu", 
      ylab = "Number of sample with otu with that proportion"
    )
  )
})

output$filter_summary_tax_before_filter = renderDataTable({
  req(get_phyloseq_data())
  return(
    summary(tax_table(get_phyloseq_data()))
  )  
})

output$filter_summary_sample_before_filter = renderDataTable({
  req(get_phyloseq_data())  
  return(
    summary(as.data.frame(sample_data(get_phyloseq_data())@.Data, col.names = sample_variables(get_phyloseq_data())))
  )
})

output$filter_plot_min_sum_after_filter = renderPlot({
  req(get_phyloseq_data())
  get_otu = function(phyobj){
    otu = otu_table(phyobj)@.Data
    tax = tax_table(phyobj)@.Data     
    if(all(rownames(tax) == rownames(otu))){
      return(otu)
    } else {
      return(t(otu))
    }     
  }
  otu = t(get_otu(physeq()))
  return(
    qplot(rowSums(otu), geom = "histogram", 
      xlab = "The amount of otu count each sample has", 
      ylab = "The number of sample with that amount of count"
    )
  )
})

output$filter_plot_zerofraction_after_filter = renderPlot({
  req(get_phyloseq_data())
  get_otu = function(phyobj){
    otu = otu_table(phyobj)@.Data
    tax = tax_table(phyobj)@.Data     
    if(all(rownames(tax) == rownames(otu))){
      return(otu)
    } else {
      return(t(otu))
    }     
  }
  otu = t(get_otu(physeq()))
  otu_zerofraction = apply(otu, MARGIN = 2, function(x) return((sum(x == 0) / length(x))))       
  
  return(
    qplot(otu_zerofraction, geom = "histogram", 
      xlab = "Proportion of zero in a otu", 
      ylab = "Number of sample with otu with that proportion"
    )
  )
})

output$filter_summary_tax_after_filter = renderDataTable({
  req(get_phyloseq_data())
  return(
    summary(tax_table(physeq()))
  )  
})

output$filter_summary_sample_after_filter = renderDataTable({
  req(get_phyloseq_data())  
  return(
    summary(as.data.frame(sample_data(physeq())@.Data, col.names = sample_variables(physeq())))
  )  
})
