# phyloseq-server-data

################################################################################
# Define the available phyloseq datasets for plotting.
################################################################################


################################################################################
# UI subset_samples expression cascade
# filter_subset_samp_expr
################################################################################
################################################################################
# The main reactive data object. Returns a phyloseq-class instance.
# This is considered the "filtered" data, used by all downstream panels,
# And generally the input to any transformation options as well
################################################################################
physeq = reactive({
  ps0 = get_phyloseq_data()
  if(input$actionb_filter == 0){
    # Don't execute filter if filter-button has never been clicked.
    if(inherits(ps0, "phyloseq")){
      return(ps0)
    } else {
      return(NULL)
    }
  }
  # Isolate all filter code so that button click is required for update
  isolate({
    if(inherits(ps0, "phyloseq")){
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
      if( !is.null(av(input$filter_samvars_selection)) ){
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
      if( input$filter_taxa_sums_threshold > 0 ){
        # OTU sums filter
        ps0 <- prune_taxa({taxa_sums(ps0) > input$filter_taxa_sums_threshold}, ps0)
      }
      if( input$filter_sample_sums_threshold > 0 ){
        # Sample sums filtering
        ps0 <- prune_samples({sample_sums(ps0) > input$filter_sample_sums_threshold}, ps0)
      }
      if(inherits(input$filter_kOverA_sample_threshold, "numeric")){
        if(input$filter_kOverA_sample_threshold > 1){
          # kOverA OTU Filtering
          flist = genefilter::filterfun(
            genefilter::kOverA(input$filter_kOverA_sample_threshold,
                               input$filter_kOverA_count_threshold, na.rm=TRUE)
          )
          koatry = try(ps0 <- filter_taxa(ps0, flist, prune=TRUE), silent = TRUE)
          if(inherits(koatry, "try-error")){
            warning("kOverA parameters resulted in an error, kOverA filtering skipped.")
          }
        }
      }
      return(ps0)
    } else {
      return(NULL)
    }
  })
})
################################################################################
################################################################################
# TRANSFORMATIONS
################################################################################
################################################################################
# Proportional transformation
################################################################################
physeqProp = reactive({
  if(is.null(output$phyloseqDataset())){
    return(NULL)
  }
  return(
    transform_sample_counts(output$phyloseqDataset(), function(x){x / sum(x)})
  )
})
################################################################################
# Regularized Log Transformation (blind, fast)
################################################################################
physeqRLog = reactive({
  if(is.null(physeq())){
    return(NULL)
  }
  ps0RLog = physeq()
  # Demo a "blind" transformation, with design formula 
  # only containing an intercept term.
  # This will often throw a warning. DESeq2 usually handles the condition fine
  dds = phyloseq_to_deseq2(ps0RLog, ~ 1)
  rld <- DESeq2::rlog(dds, blind = TRUE, fast = TRUE)
  rlogMat <- GenomicRanges::assay(rld)
  otu_table(ps0RLog) <- otu_table(rlogMat, taxa_are_rows = TRUE)
  return(ps0RLog)
})
################################################################################
# Centered Log-Ratio (CLR) Transformation
################################################################################
gm_mean = function(x, na.rm=TRUE){
  # The geometric mean, with some error-protection bits.
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}
clr = function(x, base=2){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}
physeqCLR = reactive({
  if(is.null(physeq())){
    return(NULL)
  }
  return(
    transform_sample_counts(physeq(), fun = clr)
  )
})
################################################################################
# Misc Filter-tab server code
################################################################################
# kOverA `k` Filter UI
maxSamples = reactive({
  # Create logical indicated the samples to keep, or dummy logical if nonsense input
  if(inherits(get_phyloseq_data(), "phyloseq")){
    return(nsamples(get_phyloseq_data()))
  } else {
    # Dummy response.
    return(NULL)
  }
})

# Generic Function for plotting marginal histograms
sums_hist = function(thesums=NULL, xlab="", ylab=""){
  if(is.null(thesums)){
    p = qplot(0)
  } else {
    p = ggplot(data.frame(sums=thesums), aes(x=sums))
    p = p + geom_histogram()
    p = p + xlab(xlab) + ylab(ylab) 
    p = p + scale_x_log10(labels = scales::comma)
  }
  return(p)
}
lib_size_hist = reactive({
  xlab = "Number of Reads (Counts)"
  ylab = "Number of Libraries"
  return(sums_hist(sample_sums(get_phyloseq_data()), xlab, ylab))
})
otu_sum_hist = reactive({
  xlab = "Number of Reads (Counts)"
  ylab = "Number of OTUs"
  return(sums_hist(taxa_sums(get_phyloseq_data()), xlab, ylab))    
})

################################################################################
# Component Table
################################################################################

################################################################################

# phyloseq-server-server
shinyPhyloseqServerObjectsList = ls()
rankNames = reactive({
    rankNames = as.list(rank_names(physeq(), errorIfNULL=FALSE))
    names(rankNames) <- rankNames
    return(rankNames)
  })
  variNames = reactive({
    variNames = as.list(sample_variables(physeq(), errorIfNULL=FALSE))
    names(variNames) <- variNames
    return(variNames)
  })
  vars = function(type="both", withnull=TRUE, singles=FALSE){
    if(!type %in% c("both", "taxa", "samples")){
      stop("incorrect `type` specification when accessing variables for UI.")
    }
    returnvars = NULL
    if(type=="samples"){
      if(singles){
        returnvars <- c(list(Sample="Sample"), variNames())
      } else {
        returnvars <- variNames()
      }
    }
    if(type=="taxa"){
      if(singles){
        returnvars <- c(rankNames(), list(OTU="OTU"))
      } else {
        returnvars <- rankNames()
      }
    } 
    if(type=="both"){
      # Include all variables
      if(singles){
        returnvars <- c(rankNames(), variNames(), list(OTU="OTU", Sample="Sample"))
      } else {
        returnvars <- c(rankNames(), variNames())
      }
    }
    if(withnull){
      # Put NULL first so that it is default when `select` not specified
      returnvars <- c(list("NULL"="NULL"), returnvars)
    }
    return(returnvars)
  }
  # A generic selectInput UI. Plan is to pass a reactive argument to `choices`.
  uivar = function(id, label="Variable:", choices, selected="NULL"){
    selectInput(inputId=id, label=label, choices=choices, selected=selected)
  }

# phyloseq-global
  output_phyloseq_print_html = function(physeq){
  HTML(
    paste(
      '<p class="phyloseq-print">',
      paste0(capture.output(print(physeq)), collapse=" <br/> "),
      "</p>"
    )
  )
  # Alternative tag way:
  #   do.call("p", args = c(list(class="phyloseq-print", 
  #                              sapply(c("alskfjs", "askfjls"), br, simplify = FALSE, USE.NAMES = FALSE))))
}
numericInputRow <- function(inputId, label, value, min = NA, max = NA, step = NA, class="form-control", ...){
  inputTag <- tags$input(id = inputId, type = "number", value = value, class=class, ...)
  if (!is.na(min)) 
    inputTag$attribs$min = min
  if (!is.na(max)) 
    inputTag$attribs$max = max
  if (!is.na(step)) 
    inputTag$attribs$step = step
  div(style="display:inline-block",
      tags$label(label, `for` = inputId), 
      inputTag)
}
textInputRow <- function(inputId, label, value = "", class="form-control", ...){
  div(style="display:inline-block",
      tags$label(label, `for` = inputId), 
      tags$input(id = inputId, type = "text", value = value, class=class, ...))
}
theme_blank_custom = theme_bw() + theme(
  plot.title = element_text(size = 28),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  axis.text.x      = element_blank(),
  axis.text.y      = element_blank(),
  axis.title.x     = element_blank(),
  axis.title.y     = element_blank(),
  axis.ticks       = element_blank(),
  panel.border     = element_blank()
)
shiny_phyloseq_ggtheme_list <- list(
  bl_wh = theme_bw(),
  blank = theme_blank_custom,
  thin = theme_linedraw(),
  light = theme_light(),
  minimal = theme_minimal(),
  classic = theme_classic(),
  gray = theme_gray()
)
################################################################################
# Included Data
# Define the named list of datasets to choose from
################################################################################
# Create an environment to store original loaded data
env_psdata = new.env()
# Keep server-loaded data into a special environemnt, `env_psdata`
data(list=c("GlobalPatterns", "enterotype", "esophagus"), envir = env_psdata)
load("data/exampledata.RData", envir = env_psdata)
load("data/1457_uparse.RData", envir = env_psdata)
attach(env_psdata)
# Define initial list of available datasets
datalist = list(
  closed_1457_uparse = closed_1457_uparse,
  exampledata = exampledata,
  GlobalPatterns = GlobalPatterns,
  enterotype = enterotype,
  esophagus = esophagus)
av = function(x){
  if( isTRUE(all.equal(x, "")) | isTRUE(all.equal(x, "NULL")) ){
    return(NULL)
  }
  return(x)
}
tablify_phyloseq_component = function(component, colmax=25L){
  if(inherits(component, "sample_data")){
    Table = data.frame(component)
  }
  if(inherits(component, "taxonomyTable")){
    Table = component@.Data
  }
  if(inherits(component, "otu_table")){
    if(!taxa_are_rows(component)){component <- t(component)}
    Table = component@.Data
  }
  return(Table[, 1:min(colmax, ncol(Table))])
}
# Determine available table-like components for on-screen rendering
component_options = function(physeq){
  # Initialize the return option list
  component_option_list = list("NULL"="NULL")
  # Get logical vector of components
  nonEmpty = sapply(slotNames(physeq), function(x, ps){!is.null(access(ps, x))}, ps=physeq)
  if(sum(nonEmpty)<1){return(NULL)}
  # Convert to vector of slot-name strings for non-empty components
  nonEmpty <- names(nonEmpty)[nonEmpty]
  # Cull the non-table components
  nonEmpty <- nonEmpty[!nonEmpty %in% c("phy_tree", "refseq")]
  # If no tables available, return default empty option
  if(length(nonEmpty)<1){return(component_option_list)}
  # Otherwise add to the option list and return
  compFuncString = names(phyloseq:::get.component.classes()[nonEmpty])
  if("sam_data" %in% compFuncString){
    compFuncString[compFuncString=="sam_data"] <- "sample_data"
  }
  NiceNames = c(otu_table="OTU Table",
                sample_data="Sample Data",
                tax_table = "Taxonomy Table")
  names(compFuncString) <- NiceNames[compFuncString]
  return(c(component_option_list, as.list(compFuncString)))
}
distlist <- as.list(unlist(phyloseq::distanceMethodList))
names(distlist) <- distlist
distlist <- distlist[which(!distlist %in% c("ANY"))]
