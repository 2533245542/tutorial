mistudio_miprofile = function(
	phyobj= NULL, 
	cova = NULL, conf = NULL,
	gUniFrac.alpha = NULL, gUniFrac.rarefy = NULL, pUniFrac.alpha = NULL, pUniFrac.rarefy = NULL, uwUniFrac.rarefy = NULL,
	BC.rarefy = NULL, Jaccard.rarefy = NULL,
	tree.file = NULL,
	strata = NULL,
	nperm = 1000){

	stopifnot(!is.null(phyobj))
	stopifnot(!is.null(cova)) # must select one kind of covariate (numeric or categorical)
	stopifnot(is.numeric(nperm))
	# remove NA rows in sample, but only look for NA's in selected variables
	remove_sample_na = function(
	  phyobj, 
	  cova,
	  conf,
	  strata){

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

	  find_na_index = function(selected= NULL, sample = NULL){

	    selected_split = strsplit(x = selected, split = "%mistudio_seperator%",  fixed = TRUE)[[1]]
	    stopifnot(length(selected_split) == 2)
	    value = selected_split[1] # the value users chooses, as displayed on UI
	    column = selected_split[2] # the column in which the value resides
	    if(value == column){ # only runs when user chooses the whole column
	      index = is.na(sample[, column, drop = TRUE])
	      return(index)
	    } else{
	      return(NULL)
	    }   
	  }
	  otu = get_otu(phyobj)
	  tax = get_tax(phyobj)
	  sample = get_sample(phyobj)

	  index_logical_list = lapply(c(conf, cova, strata), FUN = find_na_index, sample)
	  index_logical_list = index_logical_list[!sapply(index_logical_list, is.null)] # remove NULL lists
	  index_logical = Reduce("|", index_logical_list) # find all rows that have NA

	  sample = sample[!index_logical, ,drop = FALSE]
	  otu = otu[, !index_logical, drop = FALSE]
	  return(phyloseq(tax_table(tax), otu_table(otu, taxa_are_rows = TRUE), sample_data(sample)))
	}

	adjust_otu_dimmension = function(otu, sample){
		if(dim(otu)[1] != dim(sample)[1]) {
			return(t(otu))
		} else {
			return(otu)	
		}
	}

	get_distance = function(){

		ensure_treefile = function(){
			# if user wants to calculate tree based distances, ensure tree file is available
			if(!is.null(gUniFrac.alpha) | !is.null(pUniFrac.alpha) | !is.null(uwUniFrac.rarefy)){
				if(!is.null(phy_tree(phyobj, errorIfNULL = FALSE))){
					write.tree(phy_tree(phyobj, "temp/miprofile_Get_Tree_Dist_treefile.tre"))
					Tree.file = "temp/miprofile_Get_Tree_Dist_treefile.tre"	
				} else{
					stopifnot(!is.null(tree.file))
				}			
			}			
		}

		if(!is.null(gUniFrac.alpha) || !is.null(pUniFrac.alpha) || !is.null(uwUniFrac.rarefy)){
			ensure_treefile()
			# tree based distances included
			stopifnot(!is.null(gUniFrac.rarefy) || !is.null(pUniFrac.rarefy))
			Dist.tree = Get_Tree_Dist(OTU = otu, Tree.file = tree.file, gUniFrac.alpha = gUniFrac.alpha, 
				gUniFrac.rarefy = gUniFrac.rarefy, pUniFrac.alpha = pUniFrac.alpha, 
				pUniFrac.rarefy = pUniFrac.rarefy, uwUniFrac.rarefy = uwUniFrac.rarefy ) 
			if(is.null(BC.rarefy) && is.null(Jaccard.rarefy)){
				Dist.all = Dist.tree
			} else {
				set.seed(11)
				Dist.base = Get_Dist(OTU = otu, BC.rarefy = BC.rarefy, Jaccard.rarefy = Jaccard.rarefy)	
				Dist.all = c(Dist.base, Dist.tree)	
			}			
		} else{
			# tree based distances not included
			stopifnot(!is.null(BC.rarefy) || !is.null(Jaccard.rarefy))
			set.seed(11)
			Dist.all = Get_Dist(OTU = otu, BC.rarefy = BC.rarefy, Jaccard.rarefy = Jaccard.rarefy)	
		}
		return(Dist.all)
	}

	get_X_Z_strata = function(phyobj = NULL, cova = NULL, conf = NULL, strata = NULL){

		convert_all = function(phyobj = NULL, multiple_selected = NULL){
			
			stopifnot(!is.null(phyobj))
			stopifnot(!is.null(multiple_selected))

			convert_one = function(phyobj = NULL, selected= NULL){

				stopifnot(!is.null(phyobj))
				stopifnot(!is.null(selected))
				stopifnot(length(selected) == 1)

				# split the user selected value by sperator, assign splitted strings accordingly
				selected_split = strsplit(x = selected, split = "%mistudio_seperator%",  fixed = TRUE)[[1]]
				stopifnot(length(selected_split) == 2)

				value = selected_split[1] # the value users chooses, as displayed on UI
				column = selected_split[2] # the column in which the value resides

				sample = as.data.frame(lapply(sample_data(phyobj), as.character), stringsAsFactors = FALSE) # get sample data as a list of vectors of characters

				whole_column = !(value %in% sample[, column, drop = TRUE]) # if the whole column is selected

				if(whole_column){
					# if whole column selected 
					non_numeric = all(is.na(as.numeric(sample[, column, drop = TRUE]))) # check if selected column is actually numeric
					if(non_numeric){
						# the column is not numeric, convert to dummy var
						result = model.matrix(~., data = sample[, column, drop = FALSE])[, -1] # use [, -1] to remove intercept

					} else{
						# a numeric column
						# TODO 
						# Warning: Error in [.data.frame: undefined columns selected
						# Stack trace (innermost first):
						#     113: [.data.frame
						#     112: [
						#     111: sample_data
						#     110: [
						#     109: [
						#     108: as.matrix
						#     107: get_X_Z_strata [functions/mistudio_miprofile.R#125]
						result = as.matrix(as.numeric(sample[, column, drop = TRUE]))
					}
				} else{
					# only one sub-category selected
					result = as.matrix(ifelse(sample[column] == value, 1, 0)) # convert that value as 1, else 0
				}

				return(result)

			}

			result_list = lapply(X = multiple_selected, FUN = convert_one, phyobj = phyobj) # lapply returns a list of matrixes
			result = do.call(what = cbind, args = result_list) # cbind all the matrixes in the list

			return(result)

		}
		get_mat_keep = function(mat){
			mat_append_one = cbind(1, mat) # append one to the end of the matrix
			index_keep = qr(mat_append_one)$pivot[seq_len(qr(mat_append_one)$rank)] # get the index of columns to remove in mat_append_one
			# if all index are kept, then it means no linear dependency, thus we return the original matrix
			if(length(index_keep) == 1){ # mat is one column and all 1's
				stop("milineag/milineagc: calculating linear dependence of the cova+conf selected; all variables are linear dependent with a column of 1; please select new variables")
			}
			stopifnot(all(mat_append_one[, index_keep, drop = FALSE][, 1, drop = FALSE] == 1))
			return(mat[, (index_keep-1), drop = FALSE]) 
		}

		mat_cova = convert_all(phyobj, cova) # converts categorical to dummies, continuous to numeric columns
				
		if(is.null(conf)){
			mat_conf = NULL
		} else {
			mat_conf = get_mat_keep(convert_all(phyobj, conf))			
		}
		
		if(is.null(strata)){ #strata could only be column names
			strata = NULL
		} else {
			trimmed_name = strsplit(strata, "%mistudio_seperator%")[[1]][1] # removes the seperator to get real name			
			strata = as.matrix(sample[, trimmed_name, drop = FALSE])
		}

		result = list(X = mat_cova, Z = mat_conf, strata = strata)
		return(result)
	}

	draw_all = function(dist.all = NULL, strata = NULL){
		stopifnot(!is.null(dist.all) || !is.null(strata))

		draw_one = function(dist.one, filename, strata){
			data = as.data.frame(cmdscale(as.dist(dist.one), k = 2))
			data = cbind(data, strata =strata)
			p = ggplot(data, aes(x = V1, y = V2, col = strata)) + 
			geom_point() + 
			labs(x  = "pc1", y = "pc2")		
			ggtitle(paste(filename, "PCA Plot")) +   
			theme(plot.title = element_text(hjust = 0.5))

			ggsave(filename = paste0("temp/",filename, ".pdf"), plot = p, device = 'pdf')
			return(p)
		}
		filenames = names(dist.all)

		graph = mapply(FUN = draw_one, dist.all, filenames, MoreArgs = list(strata  = strata))
	}

	#main
	phyobj = remove_sample_na(phyobj, cova, conf, strata)
	otu = otu_table(phyobj)
	sample = sample_data(phyobj)
	otu = adjust_otu_dimmension(otu, sample)@.Data
	Dist.all = get_distance()
	X_Z_strata = get_X_Z_strata(sample, cova, conf, strata)
	X = X_Z_strata$X
	Z = X_Z_strata$Z
	strata = X_Z_strata$strata
	if(!is.null(strata)){
		strata = as.matrix(as.integer(strata), ncol = 1) # TODO delete this if josh is able to make this conversion internally
	}
	# graph = draw_all(dist.all = Dist.all, strata = as.vector(strata))
	# pval = PERMANOVA_S(Dist=Dist.all, X=X, Z=Z, strata=strata, nperm=nperm)
	# return(list(pval = pval, graph = graph))
	return(PERMANOVA_S(Dist=Dist.all, X=X, Z=Z, strata=strata, nperm=nperm))

}