mistudio_filter_otu_zerofraction = function(phyobj, fraction){
	# fracion is frmo 0 to 1 showing the otu with a 0 percentage /100 higher than that fraction will be removed
	get_otu_filter = function(phyobj, fraction){
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

		criteria_remove = apply(otu, MARGIN = 2, function(x) return((sum(x == 0) / length(x)) > fraction))		
		# find the proportion of zero in each otu, and return true if that proportion is greater than fraction
		return(otu[, !criteria_remove, drop = FALSE])		
	}
	
	get_tax_filter = function(phyobj, fraction){
		get_otu = function(phyobj){
					otu = otu_table(phyobj)@.Data
					tax = tax_table(phyobj)@.Data 		
					if(all(rownames(tax) == rownames(otu))){
						return(otu)
					} else {
						return(t(otu))
					}			
				}
		tax = tax_table(phyobj)@.Data
		otu = t(get_otu(phyobj)) 
		criteria_remove = apply(otu, MARGIN = 2, function(x) return((sum(x == 0) / length(x)) > fraction))		
		# find the proportion of zero in each otu, and return true if that proportion is greater than fraction
		return(tax[!criteria_remove, , drop = FALSE])
	}

	make_phyobj_filter = function(otu_filter,tax_filter, phyobj){
		
		return(
			phyloseq(
				otu_table(otu_filter, taxa_are_rows = FALSE),
				tax_table(tax_filter),
				sample_data(phyobj)
			)
		)
		
	}

	otu_filter = get_otu_filter(phyobj, fraction)
	browser()
	tax_filter = get_tax_filter(phyobj, fraction)
	return(make_phyobj_filter(otu_filter = otu_filter, tax_filter = tax_filter, phyobj = phyobj))


}