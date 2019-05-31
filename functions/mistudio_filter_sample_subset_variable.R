mistudio_filter_sample_subset_variable = function(phyobj, variable){
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
	get_filtered_phyobj = function(otu, tax, sample, variable){
		myfun = function(variable_one, sample){
			sample_character = as.character(sample)
			value = strsplit(variable_one, "%mistudio_seperator%")[[1]][1]
			column = strsplit(variable_one, "%mistudio_seperator%")[[1]][2]
			return(value == sample[[column]])
		}
		index_logical_list = lapply(variable, FUN = myfun, sample)
		index_logical = Reduce("|", index_logical_list)
		sample = sample[index_logical, ,drop = FALSE]
		otu = otu[, index_logical, drop = FALSE]
		return(phyloseq(tax_table(tax), otu_table(otu, taxa_are_rows = TRUE), sample_data(sample)))

	}
	otu = get_otu(phyobj)
	tax = get_tax(phyobj)
	sample = get_sample(phyobj)
	return(get_filtered_phyobj(otu, tax, sample, variable))

}