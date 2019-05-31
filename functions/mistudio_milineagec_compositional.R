# input$milineagec_categorical "LinkerPrimerSequence"
# input$milineagec_subcategorical "CCGTCAATTCMTTTRAGT%mistudio_seperator%LinkerPrimerSequence"
# input$milineagec_stratify "NECROSIS_PERCENT" "COUNTRY"
# input$milineagect_substratify "6%mistudio_seperator%NECROSIS_PERCENT"  "5%mistudio_seperator%NECROSIS_PERCENT" "1%mistudio_seperator%NECROSIS_PERCENT"  
# input$milineagec_continuous "X.SampleID"
# xlab "this is x lab"
# ylab "this is y lab"
# palette "Set1"
# theme "theme_linedraw"
# hidex FALSE
mistudio_milineagec_compositional = function(phyobj, arrange, lineage, 
	stratify, substratify, continuous,
	categorical, subcategorical,
	proportion,
	title, xlab, ylab, palette, theme, hidex
	){
	# lineag "Bacteria"
	# stratify "Age"
	# continuous "Age"
	# categorical "Age"
	make_general_plot_data = function(phyobj, lineage, stratify = NULL, continuous = NULL, categorical = NULL){		
		stopifnot(!is.null(phyobj))
		# stopifnot(length(lineage) == 1)
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
		
		# select the specified lineage	
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
				# aggregate the count for each species
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

		# convert otu count to relative abundance
		get_rel_abundance = function(subset){
			stopifnot(!is.null(subset))

			rel_helper = function(df){ # add a new column "rel" to the df

				total = sum(df$count)
				# a fail-safe in case total could be 0, which produces NaN for rel
				if(total == 0){
					df$rel = 0
				} else {
					df$rel = df$count / total	
				}				
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
	# user can stratify one plot to multiple ones
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
	  	data$rel_mean = 0
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
	filter_lineage_proportion = function(data, proportion, categorical){
		data_dt = as.data.table(data)
		if(is.null(categorical)){ # if NULL, use rel; use rel_mean otherwise
			#  rel is from 1~0 showing the rel of each lineage of each sample
			filter_tax = data_dt[, .(rel = sum(rel)), by = name][rel > proportion, 1]
			filter_tax = as.data.frame(filter_tax, stringsAsFactors = FALSE)		
			data_filtered = merge(filter_tax, data,  by.x = "name", by.y = "name", all.x = TRUE)
		} else {
			if(all(data$rel_mean == 0)){ # no subcategorical variabels selected; do nothing
				data_filtered = data
			} else {
				# rel_mean is actually count of each lineage of each sample; so we need to sum by each lineage group; then calculte the proportion
				filter_tax = data_dt[, .(rel_mean = sum(rel_mean)), by = name][, .(name, prop = rel_mean/sum(rel_mean))][prop > proportion, 1]
				filter_tax = as.data.frame(filter_tax, stringsAsFactors = FALSE)
				data_filtered = merge(filter_tax, data,  by.x = "name", by.y = "name", all.x = TRUE)	
			}			
		}
		return(data_filtered)
	}
	make_base_plot = function(data, lineage, arrange, stratify, continuous, categorical, title, xlab, ylab){								
		if(is.null(title)){
			title = arrange
		}		
		if(arrange == "categorical"){			
			if(is.null(xlab)){
				xlab = categorical
			}
			if(is.null(ylab)){
				ylab = "relative abundance"
			}
			if(is.null(stratify)){ # no stratification
				if(is.null(categorical)){ # no sampling grouping
					stopifnot(!is.null(lineage))
					if(length(lineage) == 1){
						return(
							ggplot(data, aes(x = sample_name, weight = rel, fill = name)) + 
							geom_bar(position = "fill", width = 1, color = "white") + 
							labs(title = title, x = xlab, y = ylab)						
						)
					} else{
						# more than one lineages are selected for plotting, create one plot for each 
						return(ggplot(data, aes(x = sample_name, weight = rel, fill = name)) + 
							geom_bar(position = "fill", width = 1, color = "white") + 
							facet_grid(rows = vars(parent), scales = "free") + 
							labs(title = title, x = xlab, y = ylab)						
							)
					}
					

				} else { # has sample grouping
					data[[categorical]] = as.factor(data[[categorical]]) # turn categorical column to factor, make sure bar plot x axis is discrete.
					if(length(lineage) == 1){
						return(
							ggplot(data, aes(x = eval(as.symbol(categorical)), weight = rel_mean, fill = name)) + 
							geom_bar(position = "fill", width = 1, color = "white") + 
							labs(title = title, x = xlab, y = ylab)	
						)	
					} else{
						return(
							ggplot(data, aes(x = eval(as.symbol(categorical)), weight = rel_mean, fill = name)) + 
							geom_bar(position = "fill", width = 1, color = "white") + 
							facet_grid(rows = vars(parent), scales = "free") + 
							labs(title = title, x = xlab, y = ylab)	
						)
					}
					
				}
			}
			else if(length(stratify) == 1){ # one degree of stratification
				if(is.null(categorical)){ # no sample grouping
					if(length(lineage) == 1){
						return(
							ggplot(data, aes(x = sample_name, weight = rel, fill = name)) +
								# geom_bar(position = "fill") + 
								geom_bar(position = "fill", width = 1, color = "white") +
								facet_grid(cols = vars(eval(as.symbol(stratify))), scales = "free") +
								labs(title = title, x = xlab, y = ylab)	
						)	
					} else{
						return(
							ggplot(data, aes(x = sample_name, weight = rel, fill = name)) +
								# geom_bar(position = "fill") + 
								geom_bar(position = "fill", width = 1, color = "white") +								
								facet_grid(as.formula(paste("parent", "~", stratify)), scales = "free") +
								labs(title = title, x = xlab, y = ylab)	
						)	
					}
					
					
				} else { # has sample grouping
					data[[categorical]] = as.factor(data[[categorical]]) # turn categorical column to factor, make sure bar plot x axis is discrete.
					if(length(lineage) == 1){
						return(
							ggplot(data, aes(x = eval(as.symbol(categorical)), weight = rel_mean, fill = name)) +								
								geom_bar(position = "fill", width = 1, color = "white") +
								facet_grid(cols = vars(eval(as.symbol(stratify))), scales = "free") + 
								labs(title = title, x = xlab, y = ylab)	
							)	
					} else{
						return(
							ggplot(data, aes(x = eval(as.symbol(categorical)), weight = rel_mean, fill = name)) +								
								geom_bar(position = "fill", width = 1, color = "white") +
								facet_grid(as.formula(paste("parent", "~", stratify)), scales = "free") +
								labs(title = title, x = xlab, y = ylab)	
							)
					}
					
					
				}
			} else if(length(stratify == 2)){ # two degrees of stratification
				if(is.null(categorical)){ # no sample grouping					
					if(length(lineage) == 1){
						return(
							ggplot(data, aes(x = sample_name, weight = rel, fill = name)) +							
								# geom_bar(position = "fill", width = 1, color = "white") +
								geom_bar(position = "fill") +
								facet_grid(cols = vars(eval(as.symbol(stratify[1]))), rows = vars(eval(as.symbol(stratify[2]))), scales = "free") + 
								labs(title = title, x = xlab, y = ylab)	
						)	
					} else{
						return(
							ggplot(data, aes(x = sample_name, weight = rel, fill = name)) +							
								# geom_bar(position = "fill", width = 1, color = "white") +
								geom_bar(position = "fill") +								
								facet_grid(as.formula(paste("parent", "~", stratify[1], "+", stratify[2])), scales = "free") + 
								labs(title = title, x = xlab, y = ylab)	
						)
					}
					
					
				} else { # has sample grouping
					data[[categorical]] = as.factor(data[[categorical]]) # turn categorical column to factor, make sure bar plot x axis is discrete.
					if(length(lineage) == 1){
						return(
							ggplot(data, aes(x = eval(as.symbol(categorical)), weight = rel_mean, fill = name)) +
								# geom_bar(position = "fill") + 
								geom_bar(position = "fill", width = 1, color = "white") +
								facet_grid(cols = vars(eval(as.symbol(stratify[1]))), rows = vars(eval(as.symbol(stratify[2]))), scales = "free") + 
								labs(title = title, x = xlab, y = ylab)	
							)	
					} else{
						return(
							ggplot(data, aes(x = eval(as.symbol(categorical)), weight = rel_mean, fill = name)) +
								# geom_bar(position = "fill") + 
								geom_bar(position = "fill", width = 1, color = "white") +
								facet_grid(as.formula(paste("parent", "~", stratify[1], "+", stratify[2])), scales = "free") + 
								labs(title = title, x = xlab, y = ylab)	
							)	
					}
					
					
				}
			} else {
				warning("more than one stratify selected")
				return(NULL)
			}

		} else if(arrange == "continuous"){
			if(is.null(xlab)){
				xlab = continuous
			}
			if(is.null(ylab)){
				ylab = "relative abundance"
			}
			if(is.null(stratify)){
				if(length(lineage) == 1){
					return(
						ggplot(data, aes(x = eval(as.symbol(continuous)), weight = rel, fill = name)) +
							geom_histogram(binwidth = 1, position = "fill") + 
							labs(title = title, x = xlab, y = ylab)	
						)	
				} else{
					return(
						ggplot(data, aes(x = eval(as.symbol(continuous)), weight = rel, fill = name)) +
							facet_grid(rows = vars(parent), scales = "free") + 
							geom_histogram(binwidth = 1, position = "fill") + 
							labs(title = title, x = xlab, y = ylab)	
						)	
				}
				
			}
			if(length(stratify) == 1){
				if(length(lineage) == 1){
					return(
						ggplot(data, aes(x = eval(as.symbol(continuous)), weight = rel, fill = name)) +
							geom_histogram(binwidth = 1, position = "fill") + 
							facet_grid(cols = vars(eval(as.symbol(stratify[1]))), scales = "free") + 
							labs(title = title, x = xlab, y = ylab)	
						)	
				} else{
					return(
						ggplot(data, aes(x = eval(as.symbol(continuous)), weight = rel, fill = name)) +
							geom_histogram(binwidth = 1, position = "fill") + 
							facet_grid(as.formula(paste("parent", "~", stratify)), scales = "free") +
							labs(title = title, x = xlab, y = ylab)	
						)	
				}
				
				
			} else if(length(stratify) == 2){
				if(length(lineage) == 1){
					return(
						ggplot(data, aes(x = eval(as.symbol(continuous)), weight = rel, fill = name)) +
							geom_histogram(binwidth = 1, position = "fill") + 
							facet_grid(cols = vars(eval(as.symbol(stratify[1]))), rows = vars(eval(as.symbol(stratify[2]))), scales = "free") + 
							labs(title = title, x = xlab, y = ylab)	
						)
				} else{
					return(
						ggplot(data, aes(x = eval(as.symbol(continuous)), weight = rel, fill = name)) +
							geom_histogram(binwidth = 1, position = "fill") + 
							facet_grid(as.formula(paste("parent", "~", stratify[1], "+", stratify[2])), scales = "free") + 
							labs(title = title, x = xlab, y = ylab)	
						)
				}
			} else {
				warning("more than one stratify selected")
				return(NULL)
			}
		} else {
			warning("arrange is neither categorial nor continuous")
			return(NULL)
		}
	}

	add_customized_element = function(p, arrange, palette, theme, hidex){		
		# steam line procedurally to add elements
		if(!is.null(palette)){			
			if(palette == "default"){
				p = p
			} else {
				if(palette == "1" || palette == "2" || palette == "3" || palette == "6" || palette == "7" || palette == "13"){
					palette = as.numeric(palette)
				}
				p = p + scale_fill_brewer(palette = palette)	
			}			
		}                          
		if(!is.null(theme)){
			if(theme == "default"){
				p = p
			} else if(theme == "theme_linedraw"){
				p = p + theme_linedraw()
			} else if (theme == "theme_light"){
				p = p + theme_light()
			} else if (theme == "theme_minimal"){
				p = p + theme_minimal()
			} else if (theme == "theme_classic"){
				p = p + theme_classic()
			} else if (theme == "theme_gray"){
				p = p + theme_gray()
			} else {
				warning("invalid milineaegc theme selected")
			}
		}
		if(hidex){
			p = p + theme(
				  	axis.text.x = element_blank(),				  	
				  	axis.ticks = element_blank()
				  )
			# if you want to adjust x ticks or y ticks more, check out this article:
			# http://www.sthda.com/english/wiki/ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels

		}
		# TODO: allow user to change legend as well
		# TODO: maybe allow user to select where to put legend
		p = p + labs(fill="lineage") + 
			theme(legend.position="bottom", legend.title=element_blank())
        
		return(p)

	}
########################################main########################################	

	make_data_list = function(lineage, phyobj, arrange, 
		stratify, substratify, continuous,
		categorical, subcategorical,
		proportion){
	  stopifnot(length(lineage) == 1)
		data = make_general_plot_data(phyobj, lineage, stratify, continuous , categorical) %>%
			filter_by_stratify(., stratify, substratify) %>%  
			average_by_categorical(., stratify, subcategorical, arrange) %>% 
			filter_lineage_proportion(., proportion, categorical) 
		data$parent = lineage
		return(data)
	}	
	data_list = lapply(lineage, make_data_list, phyobj, arrange, 
		stratify, substratify, continuous,
		categorical, subcategorical,
		proportion)
	data = do.call(what = rbind, args = data_list)
	p <- 
	 data %>% 
		make_base_plot(., lineage, arrange, stratify, continuous, categorical, title, xlab, ylab) %>% 
		add_customized_element(., arrange, palette, theme, hidex )
	# p <- 
	#   phyobj %>% 
	# 	make_general_plot_data(., lineage, stratify, continuous , categorical) %>%
	# 	filter_by_stratify(., stratify, substratify) %>%  
	# 	average_by_categorical(., stratify, subcategorical, arrange) %>% 
	# 	filter_lineage_proportion(., proportion, categorical) %>% 
	# 	make_base_plot(., arrange, stratify, continuous, categorical, title, xlab, ylab) %>% 
	# 	add_customized_element(., arrange, palette, theme, hidex )
	return(p)
	
}
# warning: the test below is outdated but still could used as referencial. the input "lineage" is 
# now allowed to have length greater than 1
########################################test########################################	
# # 1.1   no stratify, no categorical, no continuous
# # 1.2.1 no stratify, has categorical(all categorical), no continuous
# # 1.2.2 no stratify, has categorical(sub categorical), no continuous
# # 2.1   1stratify, no categorical, no continuous
# # 2.2.1 1stratify(all stratify), has categorical, no continuous
# # 2.2.2 1stratify(sub stratify), has categorical, no continuous
# # 3.1   2stratify, no categorical, no continuous
# # 3.2   2stratify, has categorical, no continuous
# # 3.3   2stratify(substratify), no categorical, no continuous

# # 4 no stratify, no categorical, has continuous
# # 5 1stratify, no categorical, has continuous
# # 6 2stratify, no categorical, has continuous

# # option = 1.1
# # option = 0

# test_function = function(option){
# 	# using milineagec_compositional_test.RData
# 	# continuous: Age
# 	# categorical: Sex
# 	# stratify: Race, Color

# 	load("milineagec_compositional_test.RData")
# 	phyobj = milineagec_compositional_test
# 	lineage = "Crenarchaeota"
# 	categorical = "Sex"
# 	subcategorical = c("1%mistudio_seperator%Sex", "2%mistudio_seperator%Sex", "3%mistudio_seperator%Sex")
# 	stratify = c("Race", "Color")
# 	substratify = c("1%mistudio_seperator%Race", "2%mistudio_seperator%Race", "3%mistudio_seperator%Race",
# 		"1%mistudio_seperator%Color", "2%mistudio_seperator%Color")
# 	continuous = "Age"
# 	xlab = "x lab is working"
# 	ylab = "y lab is working"
# 	palette = "Set2"
# 	theme = "theme_linedraw"
# 	hidex = FALSE
# 	if(option == "1.1" || option == "0"){

# 		mistudio_milineage_compositional(phyobj = phyobj, arrange = "categorical", lineage =lineage, 
# 			stratify = NULL, substratify = NULL, continuous = NULL,
# 			categorical = NULL, subcategorical = NULL,
# 			xlab = xlab, ylab=ylab, palette= palette, theme=theme, hidex= hidex
# 			)
# 		print("test 1.1 passed")				
# 	}
# 	if(option == "1.2.1" || option == "0"){
# 		mistudio_milineage_compositional(phyobj = phyobj, arrange = "categorical", lineage=lineage, 
# 			stratify=NULL, substratify=NULL, continuous=NULL,
# 			categorical=categorical, subcategorical=c("1%mistudio_seperator%Sex", "2%mistudio_seperator%Sex", "3%mistudio_seperator%Sex", "4%mistudio_seperator%Sex"),
# 			xlab = xlab, ylab=ylab, palette= palette, theme=theme, hidex= hidex
# 			)
# 		print("test 1.2.1 passed")						
# 	}
# 	if(option == "1.2.2" || option == "0"){
# 		mistudio_milineage_compositional(phyobj = phyobj, arrange = "categorical", lineage=lineage, 
# 			stratify=NULL, substratify=NULL, continuous=NULL,
# 			categorical=categorical, subcategorical=subcategorical,
# 			xlab = xlab, ylab=ylab, palette= palette, theme=theme, hidex= hidex
# 			)
# 		print("test 1.2.2 passed")											
# 	}
# 	if(option == "2.1" || option == "0"){
# 		mistudio_milineage_compositional(phyobj = phyobj, arrange = "categorical", lineage=lineage, 
# 			stratify=stratify[1], substratify=c("1%mistudio_seperator%Race", "2%mistudio_seperator%Race", "3%mistudio_seperator%Race","4%mistudio_seperator%Race","5%mistudio_seperator%Race"), continuous=NULL,
# 			categorical=NULL, subcategorical=NULL,
# 			xlab = xlab, ylab=ylab, palette= palette, theme=theme, hidex= hidex
# 			)
# 		print("test 2.1 passed")
# 	}
# 	if(option == "2.2.1" || option == "0"){
# 		mistudio_milineage_compositional(phyobj = phyobj, arrange = "categorical", lineage=lineage, 
# 			stratify=stratify[1], substratify=c("1%mistudio_seperator%Race", "2%mistudio_seperator%Race", "3%mistudio_seperator%Race","4%mistudio_seperator%Race","5%mistudio_seperator%Race"), continuous=NULL,
# 			categorical=categorical, subcategorical=c("1%mistudio_seperator%Sex", "2%mistudio_seperator%Sex", "3%mistudio_seperator%Sex", "4%mistudio_seperator%Sex"),
# 			xlab = xlab, ylab=ylab, palette= palette, theme=theme, hidex= hidex
# 			)
# 		print("test 2.2.1 passed")
# 	}
# 	if(option == "2.2.2" || option == "0"){
# 		mistudio_milineage_compositional(phyobj = phyobj, arrange = "categorical", lineage=lineage, 
# 			stratify=stratify[1], substratify=c("1%mistudio_seperator%Race", "2%mistudio_seperator%Race", "3%mistudio_seperator%Race"), continuous=NULL,
# 			categorical=categorical, subcategorical=c("1%mistudio_seperator%Sex", "2%mistudio_seperator%Sex", "3%mistudio_seperator%Sex", "4%mistudio_seperator%Sex"),
# 			xlab = xlab, ylab=ylab, palette= palette, theme=theme, hidex= hidex
# 			)
# 		print("test 2.2.2 passed")
# 	}
# 	if(option == "3.1" || option == "0"){
# 		mistudio_milineage_compositional(phyobj = phyobj, arrange = "categorical", lineage=lineage, 
# 			stratify=stratify, substratify=c("1%mistudio_seperator%Race", "2%mistudio_seperator%Race", "3%mistudio_seperator%Race","4%mistudio_seperator%Race","5%mistudio_seperator%Race", "1%mistudio_seperator%Color","2%mistudio_seperator%Color","3%mistudio_seperator%Color","4%mistudio_seperator%Color"), continuous=NULL,
# 			categorical=NULL, subcategorical=NULL,
# 			xlab = xlab, ylab=ylab, palette= palette, theme=theme, hidex= hidex
# 			)
# 		print("test 3.1 passed")
# 	}
# 	if(option == "3.2" || option == "0"){
# 		mistudio_milineage_compositional(phyobj = phyobj, arrange = "categorical", lineage=lineage, 
# 			stratify=stratify, substratify=c("1%mistudio_seperator%Race", "2%mistudio_seperator%Race", "3%mistudio_seperator%Race","4%mistudio_seperator%Race","5%mistudio_seperator%Race", "1%mistudio_seperator%Color","2%mistudio_seperator%Color","3%mistudio_seperator%Color","4%mistudio_seperator%Color"), continuous=NULL,
# 			categorical=categorical, subcategorical=c("1%mistudio_seperator%Sex", "2%mistudio_seperator%Sex", "3%mistudio_seperator%Sex", "4%mistudio_seperator%Sex"),
# 			xlab = xlab, ylab=ylab, palette= palette, theme=theme, hidex= hidex
# 			)
# 		print("test 3.2 passed")
# 	}
# 	if(option == "3.3" || option == "0"){
# 		mistudio_milineage_compositional(phyobj = phyobj, arrange = "categorical", lineage=lineage, 
# 			stratify=stratify, substratify=substratify, continuous=NULL,
# 			categorical=categorical, subcategorical=c("1%mistudio_seperator%Sex", "2%mistudio_seperator%Sex", "3%mistudio_seperator%Sex", "4%mistudio_seperator%Sex"),
# 			xlab = xlab, ylab=ylab, palette= palette, theme=theme, hidex= hidex
# 			)
# 		print("test 3.2 passed")
# 	}
# 	if(option == "4" || option == "0"){
# 		mistudio_milineage_compositional(phyobj = phyobj, arrange = "continuous", lineage =lineage, 
# 			stratify = NULL, substratify = NULL, continuous = continuous,
# 			categorical = NULL, subcategorical = NULL,
# 			xlab = xlab, ylab=ylab, palette= palette, theme=theme, hidex= hidex
# 			)
# 		print("test 4 passed")
# 	}
# 	if(option == "5" || option == "0"){
# 		mistudio_milineage_compositional(phyobj = phyobj, arrange = "continuous", lineage =lineage, 
# 			stratify = stratify[1], substratify = c("1%mistudio_seperator%Race", "2%mistudio_seperator%Race", "3%mistudio_seperator%Race","4%mistudio_seperator%Race","5%mistudio_seperator%Race"), continuous = continuous,
# 			categorical = NULL, subcategorical = NULL,
# 			xlab = xlab, ylab=ylab, palette= palette, theme=theme, hidex= hidex
# 			)
# 		print("test 5 passed")
# 	}
# 	if(option == "6" || option == "0"){
		
# 		mistudio_milineage_compositional(phyobj = phyobj, arrange = "continuous", lineage =lineage, 
# 			stratify = stratify, substratify = c("1%mistudio_seperator%Race", "2%mistudio_seperator%Race", "3%mistudio_seperator%Race","4%mistudio_seperator%Race","5%mistudio_seperator%Race", "1%mistudio_seperator%Color","2%mistudio_seperator%Color","3%mistudio_seperator%Color","4%mistudio_seperator%Color"), continuous = continuous,
# 			categorical = NULL, subcategorical = NULL,
# 			xlab = xlab, ylab=ylab, palette= palette, theme=theme, hidex= hidex
# 			)
# 		print("test 6 passed")
# 	}
	
# }
# # test_function("6")
# # test_function("0")



