# milineage_cova, conf, extra should be in the format of a string vector of "value%mistudio_seperatorcolumnname"
# e.g. c("1%mistudio_seperator%columnA", "Male%mistudio_seperator%Sex", "age%mistudio_seperator%age")
# TODO the writing of the program is done, but we need to test it
mistudio_milineage = function(
  phyobj = NULL,
  milineage_function = NULL,
  milineage_cova = NULL,
  milineage_conf = NULL,
  milineage_cova_extra = NULL,
  milineage_conf_extra = NULL,
  milineage_mindepth = NULL,
  milineage_nresample = NULL,
  milineage_fdralpha = NULL,
  milineage_ZI.LB = NULL,
  milineage_testtype  = NULL
){
  
  stopifnot(!is.null(phyobj))
  stopifnot(!is.null(milineage_function)) 
  stopifnot(!is.null(milineage_cova))
  
  get_OTU = function(phyobj, milineage_function){
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
  	  if(milineage_function == "QCAT_GEE"){
  	    otu = .Rarefy(otu)  
  	  }
  	  return(otu) 
  	}
  
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
  
  OTU = get_OTU(phyobj,milineage_function) # only .Rarefy for non QCAT functions(i.e. QCAT_GEE, ZIGDM)
  X = get_X(phyobj, milineage_cova, milineage_conf)
  X.index = get_X.index(phyobj, milineage_cova, milineage_conf)
  Tax = get_Tax(phyobj)
  min.depth = milineage_mindepth
  n.perm = milineage_nresample
  fdr.alpha = milineage_fdralpha
  browser()
  if(milineage_function == "QCAT"){
    return(QCAT(OTU = OTU, X = X, X.index = X.index, Tax = Tax, min.depth = min.depth, n.perm = n.perm, fdr.alpha = fdr.alpha))
    
  } else if(milineage_function == "QCAT_GEE"){
    if(is.null(milineage_cova_extra)){ # if no zero-part covaritae given, do not provide Z or Z.index
    	if(!is.null(milineage_conf_extra)){
    		warnings("milineag QCAT_GEE: ignoring zero-part confouding since zero-part covariate is not given")
    	}
    	return(QCAT_GEE(OTU = OTU, X = X, X.index = X.index, Tax = Tax, min.depth = min.depth, n.perm = n.perm, fdr.alpha = fdr.alpha))
    } else {
    	Z = get_X(phyobj, milineage_cova_extra, milineage_conf_extra)
    	Z.index = get_X.index(phyobj, milineage_cova_extra, milineage_conf_extra)
    	return(QCAT_GEE(OTU = OTU, X = X, X.index = X.index, Z = Z, Z.index = Z.index, Tax = Tax, min.depth = min.depth, n.perm = n.perm, fdr.alpha = fdr.alpha))
    }    
  } else if(milineage_function == "ZIGDM"){
    
    ZI.LB = milineage_ZI.LB
    
    if(milineage_testtype == "mean abundance"){
      
      X4mean = get_X(phyobj, milineage_cova, milineage_conf)
      X.index = get_X.index(phyobj, milineage_cova, milineage_conf)
      return(ZIGDM(OTU = OTU, X4mean = X4mean, X4disp = NULL, X4zero = NULL, test.type = "Mean", X.index = X.index, ZI.LB = ZI.LB, Tax = Tax, min.depth = min.depth, n.perm = n.perm, fdr.alpha = fdr.alpha, details = FALSE))
      
    } else if(milineage_testtype == "dispersion level"){
      
      X4disp = get_X(phyobj, milineage_cova, milineage_conf)
      X.index = get_X.index(phyobj, milineage_cova, milineage_conf)
      return(ZIGDM(OTU = OTU, X4mean = NULL, X4disp = X4disp, X4zero = NULL, test.type = "Disp", X.index = X.index, ZI.LB = ZI.LB, Tax = Tax, min.depth = min.depth, n.perm = n.perm, fdr.alpha = fdr.alpha, details = FALSE))
      
      
    } else if (milineage_testtype == "presence-absence frequency"){
      
      X4zero = get_X(phyobj, milineage_cova, milineage_conf)
      X.index = get_X.index(phyobj, milineage_cova, milineage_conf)
      return(ZIGDM(OTU = OTU, X4mean = NULL, X4disp = NULL, X4zero = X4zero, test.type = "Zero", X.index = X.index, ZI.LB = ZI.LB, Tax = Tax, min.depth = min.depth, n.perm = n.perm, fdr.alpha = fdr.alpha, details = FALSE))
      
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
  
}

########################################test########################################
# library(phyloseq)
# library(miLineage)
# load("filterbarb.RData")
# barb = filterbarb
# # select cova numeric column: BMI, Age_FecSamp, MIND_Score
# # select cova categorical column: Simple_Diagnosis, APOE4
# # select cova subcategorical： White, Female
# 
# # select conf numeric column: BMI_Age %mistudio_seperator%
# milineage_cova = c("BMI%mistudio_seperator%BMI", "Age_FecSamp%mistudio_seperator%Age_FecSamp", "MIND_Score%mistudio_seperator%MIND_Score", "Simple_Diagnosis%mistudio_seperator%Simple_Diagnosis", "APOE4%mistudio_seperator%APOE4", "White%mistudio_seperator%Ethnicity", "Female%mistudio_seperator%Sex")
# milineage_conf = "BMI_Age%mistudio_seperator%BMI_Age"
# 
# set.seed(1)
# r = mistudio_milineage(
# 	phyobj = barb,
# 	milineage_function = "ZIGDM",
# 	milineage_cova = milineage_cova,
# 	milineage_conf = milineage_conf,
# 	milineage_cova_extra = NULL,
# 	milineage_conf_extra = NULL,
# 	milineage_mindepth = 0,
# 	milineage_nresample = 1000,
# 	milineage_fdralpha = 0.05,
# 	milineage_ZI.LB = 10,
# 	milineage_testtype  = "mean abundance"
# 	)

