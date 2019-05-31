# milineagec_cova, conf, extra should be in the format of a string vector of "value%mistudio_seperatorcolumnname"
# e.g. c("1%mistudio_seperator%columnA", "Male%mistudio_seperator%Sex", "age%mistudio_seperator%age")
mistudio_milineagec = function(
	phyobj = NULL,
	milineagec_function = NULL,
	milineagec_cova = NULL,
	milineagec_conf = NULL,
	milineagec_cova_extra = NULL,
	milineagec_conf_extra = NULL,
	milineagec_mindepth = NULL,
	milineagec_nresample = NULL,
	milineagec_fdralpha = NULL,
	milineagec_id = NULL,
	milineagec_timebase = NULL,
	milineagec_permtype = NULL,
	milineagec_testtype  = NULL
	){

	stopifnot(!is.null(phyobj))
	stopifnot(!is.null(milineagec_function)) 
	stopifnot(!is.null(milineagec_cova))

	get_ID = function(phyobj, milineagec_id){
		if(is.null(milineagec_id)){
			return(NULL)
		}
		sample = as.data.frame(sample_data(phyobj)@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
			col.names = sample_data(phyobj)@names, row.names = sample_data(phyobj)@row.names) # manually reassign column and row names
		# trimmed_name = strsplit(milineagec_id, "%mistudio_seperator%")[[1]][2] # removes the seperator to get real name
		id = unlist(sample[, milineagec_id, drop = FALSE])
		return(id)
	}

	get_OTU = function(phyobj, milineagec_function){
		# unlike other get_otu function, this one is to generate OTU input for milineage so we must transpose it in the end to make row as sample and column as otu
	  otu = otu_table(phyobj)@.Data
	  tax = tax_table(phyobj)@.Data     
	  .Rarefy = function (otu.tab) 
	  {
	    # set.seed(1)
	    depth = min(rowSums(otu.tab))
	    otu.tab <- as.matrix(otu.tab)
	    ind <- (rowSums(otu.tab) < depth)
	    sam.discard <- rownames(otu.tab)[ind]
	    otu.tab <- otu.tab[!ind, ]
	    rarefy <- function(x, depth) {
	      y <- sample(rep(1:length(x), x), depth)
	      y.tab <- table(y)
	      z <- numeric(length(x))
	      z[as.numeric(names(y.tab))] <- y.tab
	      z
	    }
	    otu.tab.rff <- t(apply(otu.tab, 1, rarefy, depth))
	    rownames(otu.tab.rff) <- rownames(otu.tab)
	    colnames(otu.tab.rff) <- colnames(otu.tab)
	    return(otu.tab.rff)
	  } 
	  if(all(rownames(tax) == rownames(otu))){
	    otu = otu
	  } else {
	    otu = t(otu)
	  }     
	  otu = t(otu) # the extra line for milineage and milineagec only 
	  if(milineagec_function == "QCAT_GEE.Cluster"){
	    otu = .Rarefy(otu)  
	  }
	  return(otu) 
	}
	# 规定：这个学的一定要是每个ID group都有的，选numeric的时候可以选整个column但是categorical的话只能选subcategroy
	get_OTU.base = function(phyobj, milineagec_function, milineagec_id, milineagec_timebase){
	  get_OTU = function(phyobj, milineagec_function){
	    otu = otu_table(phyobj)@.Data
	    tax = tax_table(phyobj)@.Data     
	    .Rarefy = function (otu.tab) 
	    {
	      # set.seed(1)
	      depth = min(rowSums(otu.tab))
	      otu.tab <- as.matrix(otu.tab)
	      ind <- (rowSums(otu.tab) < depth)
	      sam.discard <- rownames(otu.tab)[ind]
	      otu.tab <- otu.tab[!ind, ]
	      rarefy <- function(x, depth) {
	        y <- sample(rep(1:length(x), x), depth)
	        y.tab <- table(y)
	        z <- numeric(length(x))
	        z[as.numeric(names(y.tab))] <- y.tab
	        z
	      }
	      otu.tab.rff <- t(apply(otu.tab, 1, rarefy, depth))
	      rownames(otu.tab.rff) <- rownames(otu.tab)
	      colnames(otu.tab.rff) <- colnames(otu.tab)
	      return(otu.tab.rff)
	    } 
	    if(all(rownames(tax) == rownames(otu))){
	      otu = otu
	    } else {
	      otu = t(otu)
	    }     
	    if(milineagec_function == "QCAT_GEE.Cluster"){
	      otu = .Rarefy(otu)  
	    }
	    
	    return(otu)
	  }
	  
	  if(is.null(milineagec_timebase)){
	    return(NULL)
	  }
	  OTU = get_OTU(phyobj, milineagec_function)
	  
	  sample = as.data.frame(sample_data(phyobj)@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
	                         col.names = sample_data(phyobj)@names, row.names = sample_data(phyobj)@row.names) # manually reassign column and row names
	  value = strsplit(milineagec_timebase, "%mistudio_seperator%")[[1]][1]
	  column = strsplit(milineagec_timebase, "%mistudio_seperator%")[[1]][2] # get the column name
	  by_otu = by(t(OTU), sample[, milineagec_id, drop = TRUE], function(x) return(x)) # note that OTU is transposed here so the row dimension allighs with sample's
	  by_sample = by(sample, sample[, milineagec_id, drop = TRUE], function(x) return(x))
	  otubase_helper = function(otu, sample, value, column){
	    otu = as.matrix(otu) # the input otu is a df because of the by() step; we convert it back to matrix here
	    if(value == column){ # numeric time base selected
	      value = min(sample[, column, drop = TRUE])
	    } else { # categorical time base selected(remember that no whole categorical is permiited to be selected for time base)
	      value = value 
	    }
	    row_index = which(sample[, column, drop = TRUE] == value)[1]
	    return(
	      do.call(rbind, replicate(nrow(otu), otu[row_index, , drop = FALSE], simplify = FALSE))
	    )
	  }
	  list_of_mat = mapply(FUN = otubase_helper, by_otu, by_sample, MoreArgs = list(value = value,  column = column), SIMPLIFY = FALSE)
	  return(
	    do.call(rbind, list_of_mat) # transpose OTU to make row as otu to be acceptable by miLineage
	  )
	}
	# load("~/Desktop/biome/_DATA/2018_0420_Barb/mistudio_barb_filter/milineagec_timebase_test.RData")
	# get_OTU.base(phyobj =milineagec_timebase_test, milineagec_function = "QCAT", milineagec_id = "ID", milineagec_timebase = "Age%mistudio_seperator%Age")
	  


	get_X = function(phyobj, cova, conf){

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

				sample = as.data.frame(sample_data(phyobj)@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
					col.names = sample_data(phyobj)@names, row.names = sample_data(phyobj)@row.names) # manually reassign column and row names

				sample = as.data.frame(lapply(sample, as.character), stringsAsFactors = FALSE) # get sample data as a list of vectors of characters

				whole_column = !(value %in% sample[, column, drop = TRUE]) # if the whole column is selected

				if(whole_column){
					# if whole column selected 
					non_numeric = suppressWarnings(all(is.na(as.numeric(sample[, column, drop = TRUE])))) # check if selected column is actually numeric
					if(non_numeric){
						# the column is not numeric, convert to dummy var
						result = model.matrix(~., data = sample[, column, drop = FALSE])[, -1] # use [, -1] to remove intercept
					} else{
						# a numeric column
						result = suppressWarnings(as.matrix(as.numeric(sample[, column, drop = TRUE])))
					}
				} else{
					# only one option selected
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
		  mat = mat_cova
		} else {
		  mat_conf = convert_all(phyobj, conf)
		  mat = cbind(mat_cova, mat_conf)
		}
		return(get_mat_keep(mat))
	}

	get_X.index = function(phyobj, cova, conf){
		# for each categorical column, do model.matrix(~., data)[-1] # remove intersect
		# for each subcategorical, do one column of 1,0 value
		# that means selecting a column at once, and selecting a all subcategories of the column would give different result
		# selecting a column at once would give one less column in X than selecting all subcategories (because model.matrix would select a baseline while selecting each subcategory would not ) 
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
					non_numeric = suppressWarnings(all(is.na(as.numeric(sample[, column, drop = TRUE])))) # check if selected column is actually numeric
					if(non_numeric){
						# the column is not numeric, convert to dummy var
						result = model.matrix(~., data = sample[, column, drop = FALSE])[, -1] # use [, -1] to remove intercept

					} else{
						# a numeric column
						result = suppressWarnings(as.matrix(as.numeric(sample[, column, drop = TRUE])))
					}
				} else{
					# only one option selected
					result = as.matrix(ifelse(sample[column] == value, 1, 0)) # convert that value as 1, else 0
				}

				return(result)

			}

			result_list = lapply(X = multiple_selected, FUN = convert_one, phyobj = phyobj) # lapply returns a list of matrixes
			result = do.call(what = cbind, args = result_list) # cbind all the matrixes in the list
			return(result)			
		}
		get_index_keep = function(mat){
			mat_append_one = cbind(1, mat) # append one to the end of the matrix
			index_keep = qr(mat_append_one)$pivot[seq_len(qr(mat_append_one)$rank)] # get the index of columns to remove in mat_append_one
			# if all index are kept, then it means no linear dependency, thus we return the original matrix
			if(length(index_keep) == 1){ # mat is one column and all 1's
				stop("milineag/milineagc: calculating linear dependence of the cova+conf selected; all variables are linear dependent with a column of 1; please select new variables")
			}
			stopifnot(all(mat_append_one[, index_keep, drop = FALSE][, 1, drop = FALSE] == 1))
			return(index_keep-1) # would give index like c(0,1,3,6,8). The 0 would not be picked when doing indexing
		}


		mat_cova = convert_all(phyobj, cova) # converts categorical to dummies, continuous to numeric columns
		if(is.null(conf)){
		  mat = mat_cova
		} else {
		  mat_conf = convert_all(phyobj, conf)
		  mat = cbind(mat_cova, mat_conf)
		}
		index_original = 1:ncol(mat_cova) # the index of cova when linear dependency is not removed
		index_final = 1:length(intersect(index_original, get_index_keep(mat)))

		return(index_final)

	}

	get_Tax = function(phyobj){
		tax = tax_table(phyobj)@.Data
		colnames(tax) = paste0("Rank", 1:ncol(tax))
		return(tax)
	}


	########################################main########################################
	ID = get_ID(phyobj, milineagec_id)	
	OTU = get_OTU(phyobj, milineagec_function)
	OTU.base = get_OTU.base(phyobj, milineagec_function, milineagec_id, milineagec_timebase)
	X = get_X(phyobj, milineagec_cova, milineagec_conf)
	X.index = get_X.index(phyobj, milineagec_cova, milineagec_conf)
	Tax = get_Tax(phyobj)
	min.depth = milineagec_mindepth
	perm.type = milineagec_permtype
	n.perm = milineagec_nresample
	fdr.alpha = milineagec_fdralpha
	test = milineagec_testtype
	if(milineagec_function == "QCAT.Cluster"){
		return(QCAT.Cluster(ID = ID, OTU = OTU, OTU.base=OTU.base, X = X, X.index = X.index, Tax=Tax, min.depth=min.depth, perm.type=perm.type, n.perm=n.perm, fdr.alpha=fdr.alpha, test=test))

	} else if(milineagec_function == "QCAT_GEE.Cluster"){

		Z = get_X(phyobj, milineagec_cova_extra, milineagec_conf_extra)
		Z.index = get_X.index(phyobj, milineagec_cova_extra, milineagec_conf_extra)
		return(QCAT_GEE.Cluster(ID = ID, OTU = OTU, OTU.base=OTU.base, X = X, X.index = X.index, Z = Z, Z.index = Z.index, Tax=Tax, min.depth=min.depth, perm.type=perm.type, n.perm=n.perm, fdr.alpha=fdr.alpha, test=test))

	} else {
		return(NULL)
	}

}

# library(phyloseq)
# library(miLineage)
# load("filterbarb.RData")
# barb = filterbarb
# # select cova numeric column: BMI, Age_FecSamp, MIND_Score
# # select cova categorical column: Simple_Diagnosis, APOE4
# # select cova subcategorical： White, Female

# # select conf numeric column: BMI_Age %mistudio_seperator%
# milineage_cova = c("MIND_Score%mistudio_seperator%MIND_Score")
# milineage_conf = c("BMI_Age%mistudio_seperator%BMI_Age", "Sequencing_Run%mistudio_seperator%Sequencing_Run")

# set.seed(1)
# mistudio_milineagec(
# 	phyobj = barb,
# 	milineagec_function = "QCAT.Cluster",
# 	milineagec_cova = milineage_cova,
# 	milineagec_conf = milineage_conf,
# 	milineagec_cova_extra = NULL,
# 	milineagec_conf_extra = NULL,
# 	milineagec_mindepth = 0,
# 	milineagec_nresample = 1000,
# 	milineagec_fdralpha = 0.05,
# 	milineagec_id = "Sequencing_Run%mistudio_seperator%Sequencing_Run",
# 	milineagec_timebase = NULL,
# 	milineagec_permtype = NULL,
# 	milineagec_testtype  = "chisq"
# 	)


