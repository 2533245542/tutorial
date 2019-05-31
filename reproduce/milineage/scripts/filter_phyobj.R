 filter_phyobj = function(
 	phyobj,
 	filter_rank_selection,
 	filter_rank,
 	filter_samvars_selection,
 	filter_samvars,
 	filter_otu_aggregate,
 	filter_otu_zerofraction,
 	filter_otu_minsum){
 	# Define generic function to access/clean variables
 	# This especially converts "NULL" to NULL
 	# av = function(x){
 	#   if( isTRUE(all.equal(x, "")) | isTRUE(all.equal(x, "NULL")) ){
 	#     return(NULL)
 	#   }
 	#   return(x)
 	# }
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
 	phyobj = phyloseq(sample_data(drop_df_factor(sample_data(phyobj))), tax_table(phyobj), otu_table(phyobj))
 	# Cascading selection filters
 	if( !is.null((filter_rank_selection)) ){
 	  keepTaxa = NULL
 	  if(!is.null(tax_table(phyobj, FALSE))){
 	    if(filter_rank == "OTU"){
 	      # OTU IDs directly
 	      keepTaxa = filter_rank_selection
 	    } else {
 	      TT = as(tax_table(phyobj), "matrix")
 	      keepTaxa = TT[, filter_rank] %in% filter_rank_selection 
 	    }
 	    if(length(keepTaxa) > 1){
 	      phyobj <- prune_taxa(keepTaxa, phyobj)
 	    } else {
 	      warning("Bad subset_taxa specification. ntaxa(phyobj) one or fewer OTUs")
 	    }
 	  }
 	}
 	if(!is.null((filter_samvars_selection)) ){
 	  keepSamples = NULL
 	  if(!is.null(sample_data(phyobj, FALSE))){
 	    if(filter_samvars == "Sample"){
 	      # Samples IDs directly
 	      keepSamples = filter_samvars_selection
 	    } else {
 	      varvec = as(get_variable(phyobj, filter_samvars), "character")
 	      keepSamples = varvec %in% filter_samvars_selection 
 	    }
 	    if(length(keepSamples) > 1){
 	      phyobj <- prune_samples(keepSamples, phyobj)
 	    } else {
 	      warning("Bad subset_taxa specification. ntaxa(phyobj) one or fewer OTUs")
 	    }
 	  }
 	}

 	if(!is.null(filter_otu_aggregate)){

 	  source("functions/mistudio_filter_otu_aggregate.R")
 	  phyobj = mistudio_filter_otu_aggregate(phyobj, filter_otu_aggregate)

 	}
 	if(filter_otu_minsum > 0 ){

 	  source("functions/mistudio_filter_otu_minsum.R")
 	  phyobj = mistudio_filter_otu_minsum(phyobj, filter_otu_minsum)
 	}
 	if(filter_otu_zerofraction > 0 ){

 	  source("functions/mistudio_filter_otu_zerofraction.R")
 	  phyobj = mistudio_filter_otu_zerofraction(phyobj, filter_otu_zerofraction)

 	}
 	
 }