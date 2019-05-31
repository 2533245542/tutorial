get_otu = function(phyobj){
	otu = otu_table(phyobj)@.Data
	tax = tax_table(phyobj)@.Data 		
	# if(dim(tax)[1] == dim(otu)[1]){
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
