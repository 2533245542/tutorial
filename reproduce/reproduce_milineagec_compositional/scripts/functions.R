get_plot_data = function(phyobj, arrange, lineage, 
	stratify, substratify, continuous,
	categorical, subcategorical,
	xlab, ylab, palette, theme, hidex){
	make_general_plot_data = function(phyobj, lineage, stratify = NULL, continuous = NULL, categorical = NULL){
			stopifnot(!is.null(phyobj))
			stopifnot(length(lineage) == 1)
			get_otu = function(phyobj){
			  otu = otu_table(phyobj)@.Data
			  tax = tax_table(phyobj)@.Data     
			  if(all(rownames(tax) == rownames(otu))){
			    return(otu)
			  } else {
			    return(t(otu))
			  }     
			}


			get_tax = function(phyobj){
				tax = tax_table(phyobj)@.Data
				return(tax)
			}

			get_sample = function(phyobj){
				sample = as.data.frame(sample_data(phyobj)@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
					col.names = sample_data(phyobj)@names, row.names = sample_data(phyobj)@row.names) # manually reassign column and row names		
				return(sample)
			}
			

			process_subset = function(otu, tax, sample, lineage, stratify= NULL, continuous = NULL, categorical = NULL){
				stopifnot(!is.null(lineage))
				stopifnot(!is.null(tax))
				stopifnot(!is.null(otu))
				stopifnot(rownames(tax) == rownames(otu)) # the two should be row aligned


				get_subset_index = function(lineage, tax){
					stopifnot(!is.null(lineage))
					stopifnot(!is.null(tax))

					filter_row = as.logical(rowSums(tax == lineage, na.rm = TRUE))
					col_sum = colSums(tax == lineage,na.rm = TRUE)

					stopifnot(! which(col_sum > 0) == length(col_sum))

					if(sum(col_sum > 0) > 1){ # warn if more than one tax rank has the lineage
						#TODO make this warning seeable from the shiny app interface
						warning(paste0("More than one tax rank has the lineage: "), lineage, "\n  Using the last rank that matches")
					} 
					filter_col = (which(col_sum > 0) + 1)[length(which(col_sum > 0))] # it is possible to have two ranks having the same lineage
														     # (I dont konw why but it happens in dataset closed_1457_uparsed for "Spirochaetes")
															 # in that case, we choose the last column that matches
					result = list(x = filter_row, y = filter_col)
					return(result)
				}

				process_otu = function(otu, tax_column){

					stopifnot(ncol(tax_column) == 1)

					process_species = function(otu_column, tax_column){
						otu_column = as.data.frame(otu_column, stringsAsFactors = FALSE)
						tax_column = as.data.frame(tax_column, stringsAsFactors = FALSE)
						df = cbind(tax_column, otu_column)
						colnames(df) = c("name", "count")
						result = aggregate(count~name, df, sum)
						return(result)

					}

					sample_names = colnames(otu)
					L = vector(mode = "list", length = length(sample_names))
					for(i in 1:length(sample_names)){
						# temp = process_otu_column(otu[i], tax)

						temp = process_species(otu[, i], tax_column)
						temp$sample_name = sample_names[i]
						
						if(!is.null(continuous)){ # if user wants to sort x by continuous
							temp[[continuous]] = sample[[continuous]][i]
						}

						if(!is.null(stratify)){ # if user wants to stratify the graph
							if(length(stratify) > 2){
								warning(paste0("More than two stratification dimensions are picked. \n Using the first two: "),stratify[1:2])
								stratify = stratify[1:2]
							}
							for (s in stratify) { # it's possible to have one or two s
								temp[[s]] = sample[[s]][i]						
							}
						}

						if(!is.null(categorical)){
							temp[[categorical]] = sample[[categorical]][i]
						}
						L[[i]] = temp
					}
					result = do.call("rbind", L)
					return(result)
				}


				index = get_subset_index(lineage, tax)

				otu_subset = otu[index$x, , drop = FALSE]
				tax_column = tax[index$x, index$y, drop = FALSE]

				result = process_otu(otu_subset, tax_column)
				return(result)
			}

			get_rel_abundance = function(subset){
				stopifnot(!is.null(subset))

				rel_helper = function(df){ # add a new column "rel" to the df
					total = sum(df$count)
					df$rel = df$count / total
					return(df)
				}
				result = do.call("rbind", # by() returns a list of df
					by(subset, subset$sample_name, FUN = rel_helper, simplify = FALSE)) 


				return(result)
			}
			otu = get_otu(phyobj)
			tax = get_tax(phyobj)
			sample = get_sample(phyobj)

			processed = process_subset(otu =otu, tax = tax, sample = sample, lineage = lineage, 
				stratify = stratify, continuous = continuous, categorical = categorical)
			rel = get_rel_abundance(processed)	
			return(rel)

		}
		# substratify "6%mistudio_seperator%NECROSIS_PERCENT" "20%mistudio_seperator%NECROSIS_PERCENT"
		filter_by_stratify = function(data, stratify, substratify){
			if(is.null(stratify)){
				return(data)
			}
			if(length(stratify) == 1){
				# for each substratify, split by mistudio_seperator to get the value and the corresponding column
				# for each of these pairs, do data[column] == value to find the rows that need to be kept
				# OR all the result to get a boolean vector with the length equal to the nrow of data
				# use that boolean vector to filter data by row and return
				rows_keep_criteria = 
					Reduce('|', 
						lapply(
							strsplit(substratify, "%mistudio_seperator%"), 
							function(valuecolumn, data){
								data[[valuecolumn[2]]] == valuecolumn[1]
							}, 
							data = data)
					)
				if(all(rows_keep_criteria == FALSE)){
					warning("")
				}
				
				return(data[rows_keep_criteria, , drop = FALSE])
			}
			if(length(stratify) == 2){			
				rows_keep_criteria1 = 
					Reduce('|', 
						lapply(
							strsplit(grep(stratify[1], substratify, value = TRUE), "%mistudio_seperator%"), 
							function(valuecolumn, data){
								data[[valuecolumn[2]]] == valuecolumn[1]
							}, 
							data = data)
					)
				rows_keep_criteria2 = 
					Reduce('|', 
						lapply(
							strsplit(grep(stratify[2], substratify, value = TRUE), "%mistudio_seperator%"), 
							function(valuecolumn, data){
								data[[valuecolumn[2]]] == valuecolumn[1]
							}, 
							data = data)
					)
				rows_keep_criteria0 = rows_keep_criteria1 & rows_keep_criteria2
				return(data[rows_keep_criteria0, , drop = FALSE])			
			}
			
			
		}

	# subcategorical "6%mistudio_seperator%NECROSIS_PERCENT" "9%mistudio_seperator%NECROSIS_PERCENT"
		average_by_categorical = function(data, stratify, subcategorical, arrange){
			  if(arrange == "continuous"){
			    return(data)
			  }
			  if(is.null(subcategorical)){
			    return(data)
			  }
			  
			  categorical = strsplit(subcategorical, "%mistudio_seperator%")[[1]][2]
			  rows_keep_criteria = 
			    Reduce('|', 
			           lapply(
			             strsplit(subcategorical, "%mistudio_seperator%"), 
			             function(valuecolumn, data){
			               data[[valuecolumn[2]]] == valuecolumn[1] # it is possible to have valuecolumn[1] be "12" etc. But no need to convert data[[valuecolumn[2]]] to be character because in R, "2" == 2 is TRUE
			             }, 
			             data = data)
			    )
			  filter_data = data[rows_keep_criteria, , drop = FALSE] # filter by subcategorical
			  if(length(stratify) == 0){
			    list_of_row = by(filter_data, list(filter_data$name, filter_data[[categorical]]), function(x) {
			      a_row = x[1, , drop = FALSE] # everything in x except ID, count and rel is unique. We will modify count and rel later, and ID does not matter because it would never be used 
			      a_row$count = sum(x$count)
			      a_row$rel_mean = sum(x$count)
			      return(a_row)}, simplify = FALSE)
			    return(
			      do.call("rbind", list_of_row[!sapply(list_of_row, is.null)])  
			    )
			  } else if(length(stratify) == 1){
			    list_of_row = by(filter_data, list(filter_data$name, filter_data[[stratify[1]]],  filter_data[[categorical]]), function(x) {
			      a_row = x[1, , drop = FALSE]
			      a_row$count = sum(x$count)
			      a_row$rel_mean = sum(x$count)
			      return(a_row)}, simplify = FALSE)
			    return(
			      do.call("rbind", list_of_row[!sapply(list_of_row, is.null)])  
			    )
			  } else if(length(stratify) == 2){
			    # group by name, stratify[1], stratify[2], categorical and then aggregate the count
			    list_of_row = by(filter_data, list(filter_data$name, filter_data[[stratify[1]]], filter_data[[stratify[2]]], filter_data[[categorical]]), function(x) {
			      a_row = x[1, , drop = FALSE]
			      a_row$count = sum(x$count)
			      a_row$rel_mean = sum(x$count)
			      return(a_row)}, simplify = FALSE)
			    return(
			      do.call("rbind", list_of_row[!sapply(list_of_row, is.null)])  
			    )
			    
			  } else {
			    warning("average_by_categorical: more than two stratification selected")
			  }
			}
	general_plot_data = make_general_plot_data(phyobj, lineage, stratify, continuous , categorical)
	stratify_data = filter_by_stratify(general_plot_data, stratify, substratify)
	return(average_by_categorical(stratify_data, stratify, subcategorical, arrange))
}
 