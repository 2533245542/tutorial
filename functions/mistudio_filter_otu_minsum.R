mistudio_filter_otu_minsum = function(phyobj, minsum){
	get_otu_filter = function(phyobj, minsum){

		get_otu = function(phyobj){
					otu = otu_table(phyobj)@.Data
					tax = tax_table(phyobj)@.Data 		
					if(all(rownames(tax) == rownames(otu))){
						return(otu)
					} else {
						return(t(otu))
					}			
				}
		otu = t(get_otu(phyobj))
		criteria = apply(otu, MARGIN = 1, function(x) (sum(x) >= minsum)) # find all rows that have otu sum >= minsum		
		return(otu[criteria, , drop = FALSE])
	}

	get_sample_filter = function(phyobj, minsum){
		get_otu = function(phyobj){
					otu = otu_table(phyobj)@.Data
					tax = tax_table(phyobj)@.Data 		
					if(all(rownames(tax) == rownames(otu))){
						return(otu)
					} else {
						return(t(otu))
					}			
				}
		sample = as.data.frame(sample_data(phyobj)@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
			col.names = sample_data(phyobj)@names, row.names = sample_data(phyobj)@row.names) # manually reassign column and row names		

		otu = t(get_otu(phyobj)) 
		criteria = apply(otu, MARGIN = 1, function(x) (sum(x) >= minsum)) # find all rows that have otu sum > minsum		
		return(sample[criteria, , drop = FALSE])
	}

	make_phyobj_filter = function(otu_filter,sample_filter, phyobj){
		
		return(
			phyloseq(
				otu_table(otu_filter,taxa_are_rows = FALSE),
				tax_table(phyobj),
				sample_data(sample_filter)
			)
		)
		
	}

	otu_filter = get_otu_filter(phyobj, minsum)
	sample_filter = get_sample_filter(phyobj, minsum)
	# hist(rowSums(otu_filter))
	return(make_phyobj_filter(otu_filter = otu_filter, sample_filter = sample_filter, phyobj = phyobj))

}

# library(phyloseq)
# load("newbarb.RData")
# mistudio_filter_otu_minsum(barb, minsum = 60000)