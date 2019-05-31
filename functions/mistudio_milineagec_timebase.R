# phyobj = GlobalPatterns, ID = "Group"
mistudio_milineagec_timebase = function(phyobj, ID){
    get_sample = function(phyobj){
  	  sample = as.data.frame(sample_data(phyobj)@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
  	  	col.names = sample_data(phyobj)@names, row.names = sample_data(phyobj)@row.names) # manually reassign column and row names		
    }

	make_timebase = function(phyobj, sample, ID){

		select_common_sample = function(sample, ID){
			# group sample by ID, return the grouped data frames
	  		list_of_dfs = by(sample, sample[, ID, drop = TRUE], function(x) {return(x)}, simplify = FALSE)
	  		
	  		filter_list = lapply(1:ncol(sample), 
					  			FUN = function(x){
					  				list_of_sets = lapply(list_of_dfs, FUN = '[', x) # subset the xth column of each df
					  				if(suppressWarnings(!any(is.na(as.numeric(as.character(unlist(lapply(list_of_sets, unlist)))))))){
					  					# when the xth column is numeric, the column is acceptable for time base selection, return the whole column, 
					  					return(as.numeric(as.character(unlist(lapply(list_of_sets, unlist)))))
					  				} else {
					  					# when the xth column is categorial, do intersect to find common categories and return the common categories 
					  					return(
					  						Reduce(intersect, lapply(list_of_sets, FUN = function(x) unlist(x))) # convert each list_of_set from df to a vector, intersect them, return the overlapping elements
					  					)					  					
					  				}
					  				
				  			})
	  		names(filter_list)  = colnames(sample)
	  		filter_list[sapply(filter_list, function(x) return(length(x) == 0))] <- NULL # remove all the NULL's in filter_list
	  		filter_list[ID] = NULL # remove ID from time base options
	  		return(filter_list)	
		}

		get_filter_sample = function(phyobj, sample, filter_list){
			
			mistudio_milineage_cova = function(phyobj = NULL){

				stopifnot(!is.null(phyobj))

				filter_list = function(vector = NULL, name = NULL){
					# unique each incoming vector, return 

					stopifnot(!is.null(vector))
					stopifnot(!is.null(name))

					uniqued_vector = unique(vector)

					all_numeric = suppressWarnings(!any(is.na(as.numeric(uniqued_vector))))
								  # cast uniqued_vector to numeric, check if result is all na
								  # if so, it means uniqued_vector is all character and thus
								  # needs to be abandoned					
								  # this line is meant to have warnings so I suppress them to not confuse users

					# if(length(uniqued_vector) == 1){
					# 	# only one unique value cant do analysis
					# 	return(NULL)
					# }

					# if(length(uniqued_vector) > 10 && !all_numeric){
					# 	# more than 10 unique values and they not all numeric
					# 	return(NULL)
					# }

					if(all_numeric){
						vector = vector[0]
						attr(vector, "info") <- "numeric" # give data names that reflect numeric or categorical
						return(vector)
					}

					if(!all_numeric){
						attr(uniqued_vector, "info") <-  "categorical" # give data names that reflect numeric or categorical
						return(uniqued_vector)
					}
					
					return(uniqued_vector)

				}

				append_extra = function(list = NULL, old_names = NULL){
					# make the strings of the vector of the list a named string
					# the sole purpose of this is it is possible for two cova columns having duplicated names
					# e.g.  `East Coast` = c("NY" = "NY-East Coast"),
				    #       `West Coast` = c("NY" = "NY-West Coast")
				    # Notice that both East Coast and West Coast have NY. Also be informed that only "NY" the attribute
				    # is displayed to users, but when users select, the value is the one to be returned, i.e. "NY-East Coast".
					# By doing so, we can know which NY the user chooses even though they look exactly the same in UI display.
					stopifnot(!is.null(list))

					per_str = function(string, title) paste(string, title, sep = "%mistudio_seperator%", collapse = NULL) # cant use paste0 cuz it does not have the sep arg
					per_vec = function(vec, title) {			
						vec = append(vec, title, 0) # insert the column name to be the first element
						sapply(X = vec, FUN = per_str, title = title, simplify = TRUE, USE.NAMES = TRUE)
					}
						# important to know that simplify = TRUE to not return list, USE.NAMES = true to set correct attributes
					per_list = function(list) mapply(FUN = per_vec, list, old_names, SIMPLIFY = FALSE)
					result = per_list(list) # a list consists of vectors consists of strings
					return(result)

				}


				indent_option = function(list = NULL, name = NULL, leading_space = NULL){
					list = append(list, name, 0) # make the first element of the list the name of itself
					indented = sapply(list[2:length(list)], FUN = function(x) paste0(leading_space, x)) # indent 2nd to the last element
					combined = c(name, indented) # combine the first and the rest elements to a vector, return
					return(combined)
				}

				data = lapply(sample_data(phyobj), as.character) # get sample data as a list of vectors of characters
				filtered_data = mapply(filter_list, data, names(data),  SIMPLIFY = TRUE, USE.NAMES = TRUE) # get unique elements and apply filters and make list name reflect numeric or categorical
				
				cleaned_data = filtered_data[!sapply(filtered_data, is.null)] # remove NULL lists
				
				# old_names = names(cleaned_data)
				
				info = sapply(cleaned_data, function(x) attr(x, "info")) # give data names that reflect numeric or categorical
				
				numeric_col = names(info[info == "numeric"]) # find col names that are numeric columns
				categorical_col = names(info[info == "categorical"]) # find col names that are categorical columns

				indented_data = mapply(indent_option, cleaned_data[categorical_col], categorical_col, # make the col name the first element of each list
					MoreArgs = list(leading_space = "--"), SIMPLIFY = FALSE)          # and indent the 2nd to last elements to make them 
																					  # distinguishable from col names
				
				appended_data = append_extra(cleaned_data[categorical_col], categorical_col) # add %mistudio_seperator% 
				
				result_categorical = unlist(appended_data, use.names = FALSE) # unlist the categorical options to make it a vector
				names(result_categorical) = unlist(indented_data, use.names = FALSE) # give it names for display to user purpose
				
				result_numeric = paste0(numeric_col, "%mistudio_seperator%", numeric_col) # add mistudio_seperator, just to make it align with categorical data
				names(result_numeric) = numeric_col # for display purpose
				
				if(length(numeric_col) == 0) # make sure no numeric to display when there is no numeric column
					result_numeric = NULL
				if(length(categorical_col) == 0) # make sure no categorical to display when there is no categorical column
					result_categorical = NULL

				result = list(numeric = result_numeric, categorical = result_categorical)
				return(result)

			}

			########## start ##########
			if(length(filter_list) == 0){
				return("NULL")
			}
			filter_sample <- 
				filter_list %>%
				lapply(., FUN = function(x, len){ # extend x, so that we have a data frame as qualified input to mistudio_milineage_cova
			  					if(is.numeric(x[1])){
			  						model = numeric(len)
			  					} else {
			  						model = character(len)
			  					}
			  					if(length(x) == length(model)){
			  						model = x
			  					} else if(length(x) < length(model)){			  						
			  						model[1:length(x)] = x
			  						model[(length(x) + 1) : len] = x[1]
			  					} else {
			  						stop("milineagec_timebase creating a dataframe longer than sample")
			  					}
			  					return(model)
			  				},
	  						len = nrow(sample)) %>%
				as.data.frame(., stringsAsFactors = FALSE, col.names = names(filter_list), row.names = rownames(sample)) %>%
				phyloseq(sample_data(.), tax_table(phyobj), otu_table(phyobj)) %>%
				mistudio_milineage_cova(.)
			return(filter_sample)
		}  
		filter_list = select_common_sample(sample, ID)
		filter_sample = get_filter_sample(phyobj, sample, filter_list)
		return(filter_sample)
		
	}
	sample = get_sample(phyobj)
	return(make_timebase(phyobj, sample, ID))
}
########################################test########################################
# load("milineagec_timebase_test.RData")
# mistudio_milineagec_timebase(phyobj = milineagec_timebase_test, ID = "ID")
