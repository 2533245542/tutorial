mistudio_filter_sample_subset_subject = function(phyobj, subject){
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
	get_filtered_phyobj = function(otu, tax, sample, subject){
		sample = sample[subject, ,drop = FALSE]
		otu = otu[, subject, drop = FALSE]
		return(phyloseq(tax_table(tax), otu_table(otu, taxa_are_rows = TRUE), sample_data(sample)))
	}
	otu = get_otu(phyobj)
	tax = get_tax(phyobj)
	sample = get_sample(phyobj)
	return(get_filtered_phyobj(otu, tax, sample, subject))
}