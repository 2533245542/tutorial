mistudio_filter_otu_aggregate = function(phyobj, filter_rank){
	get_tax = function(phyobj){
		tax = tax_table(phyobj)@.Data
		return(tax)
	}

	get_otu = function(phyobj){
		otu = otu_table(phyobj)@.Data
		tax = tax_table(phyobj)@.Data 		
		if(all(rownames(tax) == rownames(otu))){
			return(otu)
		} else {
			return(t(otu))
		}			
	}


	get_aggregate = function(otu, tax, filter_rank){
		get_otu_filter = function(otu, tax, filter_rank){
			otu_filter = t(aggregate(t(otu), list(tax[, filter_rank, drop = TRUE]), FUN = sum)[-1]) # aggregate based on unique tax[filter_rank]
																			   # then sum up all rows in a group, and [-1] to
																			   # remove the group variable column
			colnames(otu_filter) = paste0(filter_rank, na.omit(unique(tax[, filter_rank, drop = TRUE])))
																				  # unique_value[n] is possible to be "", if so
																				  # for row name it will be converted to something like "4" where the row number is used
																				  # for col name it will be "V4" causing tax, otu mismatch when make phyobj
																				  # so have to add something to make sure it is not ""
			return(otu_filter)
		}
		get_tax_filter = function(otu, tax, filter_rank){
			unique_values = unique(tax[, filter_rank, drop = TRUE])
			tax_filter = matrix(character(), nrow = length(unique_values), ncol = which(colnames(tax) == filter_rank)) # create an empty matrix
			for(i in 1:length(unique_values)){ # loop over the index of length(unique values in tax[filter_rank])
				tax_filter[i, ] = tax[which(tax[, filter_rank, drop = TRUE] == unique_values[i])[1], ,drop = FALSE][1, 1:which(colnames(tax) == filter_rank), drop = FALSE] # subset the row to fit 
								  	  # the index of first row in filter_rank that is equal #######
									  # unique value to the ith ################################### 
								  # select that row ############################################################# #
																												   # subset that row from the beginning to the column #####
																												   # of filter_rank #######################################															

				# find a group of rows in tax[,filter_rank] that 
				# matches the ith  unque value, and add the first row of the group to matrix
				# at the mean time, subset the row to fit tax_filter, and add it to tax_filter
			}
			colnames(tax_filter) = colnames(tax)[1:which(colnames(tax) == filter_rank)] # select from 1 to the filter_rank's column names
			rownames(tax_filter) = paste0(filter_rank, unique_values) # unique_value[n] is possible to be "", if so
																	  # for row name it will be converted to something like "4" where the row number is used
																	  # for col name it will be "V4" causing tax, otu mismatch when make phyobj
																	  # so have to add something to make sure it is not ""
			return(tax_filter)
		}

		make_phyobj_filter = function(otu_filter, tax_filter, phyobj){
			
			return(
				phyloseq(
					otu_table(otu_filter, taxa_are_rows = FALSE),
					tax_table(tax_filter),
					sample_data(phyobj)
				)
			)
			
		}

		otu_filter = get_otu_filter(otu, tax, filter_rank)
		tax_filter = get_tax_filter(otu, tax, filter_rank)
		return(make_phyobj_filter(otu_filter = otu_filter, tax_filter = tax_filter, phyobj = phyobj))
	}
########################################main########################################
	otu = t(get_otu(phyobj))
	# a = otu
	# hist(colSums(a))
	tax = get_tax(phyobj)
	return(get_aggregate(otu = otu, tax = tax, filter_rank = filter_rank))
}
# load("newbarb.RData")
# mistudio_filter_otu_aggregate(phyobj = barb, filter_rank = "ta5")

