mistudio_filter_tax_subset = function(phyobj, rank, taxa){
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
	get_filtered_phyobj = function(otu, tax, sample, rank, taxa){
		if(!is.null(taxa)){
			stopifnot(taxa != "NA")
			tax_row_compared_list = lapply(taxa, function(x, tax)return(tax == x), tax)
			tax_row_index = as.logical(rowSums(Reduce('|', tax_row_compared_list), na.rm = TRUE))
			tax = tax[tax_row_index, , drop = FALSE]
			otu = otu[tax_row_index, , drop = FALSE]
		}
		if(!is.null(rank)){
			rownames_tax = rownames(tax)
			tax_col_list = lapply(rank, function(x, tax) return(tax[, x, drop = TRUE]), tax)
			tax = do.call(cbind, tax_col_list)
			# colnames(tax) = rank
			colnames(tax) = paste0("Rank", 1:length(rank))
			rownames(tax) = rownames_tax
		}

		return(phyloseq(tax_table(tax), otu_table(otu, taxa_are_rows = TRUE), sample_data(sample)))
	}
	otu = get_otu(phyobj)
	tax = get_tax(phyobj)
	sample = get_sample(phyobj)
	return(get_filtered_phyobj(otu, tax, sample, rank, taxa))
}