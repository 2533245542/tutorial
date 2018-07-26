mistudio_milineage_compositional = function(phyobj, lineage, sortby = NULL, stratifyby = NULL){
	stopifnot(!is.null(phyobj))
	stopifnot(length(lineage) == 1)
	if(!is.null(sortby)) stopifnot(length(sortby) == 1)
	if(!is.null(stratifyby)) stopifnot(length(stratifyby) == 1)

	adjust_otu_dimmension = function(otu, tax, sample){

		otu_dim = dim(otu)
		tax_dim = dim(tax)
		sample_dim = dim(sample)

		if(otu_dim[2] == tax_dim[1]) otu = as.data.frame(t(otu), stringsAsFactors = FALSE)

		otu_dim_new = dim(otu)

		stopifnot(otu_dim_new[1] == tax_dim[1])
		stopifnot(otu_dim_new[2] == sample_dim[1])

		return(otu)
	}

	process_subset = function(otu, tax, sample, lineage, sortby = NULL, stratifyby= NULL){

		stopifnot(!is.null(lineage))
		stopifnot(!is.null(tax))
		stopifnot(!is.null(otu))
		stopifnot(rownames(tax) == rownames(otu)) # the two should be row aligned


		get_subset_index = function(lineage, tax){
			stopifnot(!is.null(lineage))
			stopifnot(!is.null(tax))

			filter_row = as.logical(rowSums(tax == lineage, na.rm = TRUE))
			col_sum = colSums(tax == lineage,na.rm = TRUE)

			stopifnot(!sum(col_sum > 0) > 1)
			stopifnot(! which(col_sum > 0) == length(col_sum))
			filter_col = (which(col_sum > 0) + 1)
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
				if(!is.null(sortby))temp$sortby = sample[[sortby]][i]
				if(!is.null(stratifyby))temp$stratifyby = sample[[stratifyby]][i]
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

	plot_compositional= function(data, sortby= NULL, stratifyby = NULL ){
		
		if(sortby == "NULL") sortby = NULL # "" means sortby is not selected (b/c of Shiny's restriction)
		if(stratifyby == "NULL") stratifyby = NULL
		
		if(is.null(sortby) & is.null(stratifyby)){
			result = ggplot(data, aes(x = sample_name, weight = rel, fill = name)) +
			  geom_bar(position = "fill")
		} else if(!is.null(sortby) & is.null(stratifyby)){
			result = ggplot(data, aes(x = sortby, weight = rel, fill = name)) +
			  geom_histogram(binwidth = 1, position = "fill")

		} else if(is.null(sortby) & !is.null(stratifyby)){
			result = ggplot(data, aes(x = sample_name, weight = rel, fill = name)) +
			  geom_bar(position = "fill") + 
			  facet_wrap(~stratifyby, ncol = 1)

		} else { #both sortby and stratifyby are not NULL
			result = ggplot(data, aes(x = sortby, weight = rel, fill = name)) +
			  geom_histogram(binwidth = 1, position = "fill") + 
			  facet_wrap(~stratifyby, ncol = 1)
		}

		return(result)
	}

	otu = as.data.frame(otu_table(phyobj), stringsAsFactors = FALSE)
	tax = as.data.frame(tax_table(phyobj), stringsAsFactors = FALSE)
	sample = as.data.frame(sample_data(phyobj), stringsAsFactors = FALSE)

	otu = adjust_otu_dimmension(otu, tax, sample)

	processed = process_subset(otu =otu, tax = tax, sample = sample, lineage = lineage, sortby = sortby, 
		stratifyby = stratifyby)
	rel = get_rel_abundance(processed)

	plot_compositional(rel, sortby = sortby, stratifyby = stratifyby)
}