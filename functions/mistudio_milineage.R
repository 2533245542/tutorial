# library(phyloseq)
# library(caret)
# library(miLineage)

mistudio_milineage = function(phyobj= NULL, cov_n= NULL, cov_c = NULL, con_n= NULL, con_c= NULL, 
	n.resample = 1000, fdr.alpha = 0.05, min.depth = 0){

	library(phyloseq)

	stopifnot(!is.null(phyobj))
	stopifnot(!is.null(cov_n) | !is.null(cov_c)) # must select one kind of covariate (numeric or categorical)

	stopifnot(is.numeric(n.resample))
	stopifnot(is.numeric(fdr.alpha))

	# > data.real$OTU.real %>% dim()
	# [1] 96 80
	# > data.real$Tax.real %>% dim()
	# [1] 80  6
	# > data.real$covariate.real %>% dim()

	milineage = function(phyobj, cov_n, cov_c, con_n, con_c, 
		n.resample = 1000, fdr.alpha = 0.05, min.depth = 0){

		adjust_otu_dimmension = function(otu, tax, sample){

			otu_dim = dim(otu)
			tax_dim = dim(tax)
			sample_dim = dim(sample)

			if(otu_dim[1] == tax_dim[1]) otu = t(otu)

			otu_dim_new = dim(otu)
			stopifnot(otu_dim_new[2] == tax_dim[1])
			stopifnot(otu_dim_new[1] == sample_dim[1])

			return(otu)
		}

		get_sample.QCAT_GEE_index = function(sample = NULL, cov_n = NULL, cov_c = NULL, con_n = NULL, con_c = NULL){

			get_dummy = function(sample, cov_c, con_c){
				data = as.data.frame(x = sample@.Data, row.names = sample@row.names, col.names = sample@names)
				model = dummyVars(~., data = data)
				dummy = predict(model, data)
				index_cov = sapply(cov_c, function(X) grep(X, colnames(dummy)))
				index_con = sapply(con_c, function(X) grep(X, colnames(dummy)))
				
				print(as.numeric(index_cov))
				print(dummy[, as.numeric(index_cov)])

				result = list(cov = as.data.frame(dummy[, as.numeric(index_cov)]), # cast to data.frame in case index_cov
					con = as.data.frame(dummy[, as.numeric(index_con)]))		   # has length 1 and returned as a vector
				return(result)													   # That would make ncol(dummy$cov) == 0
			}

			remove_linear_independent_columns = function(m){
				result = m[, qr(m)$pivot[seq_len(qr(m)$rank)]] # StackOverflow solution
				return(result)
			}

			if(is.null(cov_c) & is.null(con_c)){ # no categorical variable selected

				if(is.null(con_n)){ # no numeric confounding column selected
					mat = as.matrix(sample[, cov_n])
				} else{
					mat = as.matrix(cbind(sample[, cov_n], sample[, con_n]))	
				}

				index = 1:length(cov_n)
				result = list(mat = mat, index = index)
				return(result)
			}
			
			dummy = get_dummy(sample, cov_c, con_c)

			# combine categorical and numeric columns
			if(is.null(cov_n)){
				if(is.null(con_n)){
					mat = cbind(dummy$cov, dummy$con)
				}
				else{
					mat = cbind(dummy$cov, dummy$con, sample[ ,con_n])
				}
			} else{
				if(is.null(con_n)){
					mat = cbind(dummy$cov, sample[ ,cov_n], dummy$con)
				}
				else{
					mat = cbind(dummy$cov, sample[ ,cov_n], dummy$con, sample[ ,con_n])
				}
				
			}

			mat = as.matrix(remove_linear_independent_columns(mat))

			index = 1:(ncol(dummy$cov) + length(cov_n))
			result  = list(mat = mat, index = index)
			return(result)
		}

		.Rarefy = function (otu.tab) 
		{
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
		

		library(caret)
		library(miLineage)


		# Get sample, otu, tax from physeq obj
		otu = otu_table(phyobj)
		tax = tax_table(phyobj)
		sample = sample_data(phyobj)

		otu = adjust_otu_dimmension(otu, tax, sample) # transpose otu if wrong dimmension

		set.seed(11) #for exampledata producibility

		otu.QCAT_GEE = .Rarefy(otu)
		tax.QCAT_GEE = tax@.Data
		colnames(tax.QCAT_GEE) = paste0("Rank", 1:ncol(tax.QCAT_GEE))
		result_get_sample.QCAT_GEE_index = get_sample.QCAT_GEE_index(sample, cov_n, cov_c, con_n, con_c)
		sample.QCAT_GEE = result_get_sample.QCAT_GEE_index$mat
		indices.cova.QCAT_GEE = result_get_sample.QCAT_GEE_index$index
		browser()

		QCAT_GEE(otu.QCAT_GEE, sample.QCAT_GEE, indices.cova.QCAT_GEE,
		         sample.QCAT_GEE, indices.cova.QCAT_GEE, tax.QCAT_GEE,
		         n.perm = n.resample, fdr.alpha=fdr.alpha, min.depth = min.depth)
	}

	graphlan = function(phyobj, result){

		if(is.null(result$sig.lineage)){
			warning("No significant lineage found in miLineage")
			return()
		}

		
		## function for write tax
		writetax = function(phyobj){
		  tax = tax_table(phyobj)
		  ## unique tax
		  tax = unique(tax)
		  fileinfo = file('guide.txt')

		  i = 1
		  while(i <= nrow(tax)){
		    # parse NA
		    names = as.character(tax[i,])
		    index = as.logical(is.na(names))
		    names = names[!index]
		    ## set string
		    writestr = paste0(names, collapse = ".")
		    
		    write(writestr, 'guide.txt', append = TRUE)
		    i = i + 1
		  }
		  close(fileinfo)
		  
		}

		writeannot = function(results){
		  
		  # ************************test************************ #
		  # results = c("Erysipelotrichaceae","Erysipelotrichaceae",
		  #             "Neisseriaceae", "Pseudomonadaceae" )
		  # ************************test************************ #
		  
		  results = unique(results)
		  if(length(results) == 0){
		    return()
		  }
		  ########### Start Ugly Work Around Method For Just exampledata########### 
		  # results = substr(results, 4, nchar(results))
		  ########### End Ugly Work Around Method For Just exampledata########### 
		  fileinfo = file('annot.txt')
		  cat(file = fileinfo)  
		  
		  write("branch_bracket_depth\t0.8", 'annot.txt', append = TRUE)
		  write("branch_bracket_width\t0.25", 'annot.txt', append = TRUE)
		  write("branch_thickness\t0.5", 'annot.txt', append = TRUE)
		  write("clade_marker_size\t1", 'annot.txt', append = TRUE)
		  write("clade_marker_edge_color\t#555555", 'annot.txt', append = TRUE)
		  write("clade_marker_edge_width\t0.1", 'annot.txt', append = TRUE)
		  write("", 'annot.txt', append = TRUE)

		  # Get rainbow colors 
		  colors.original = rainbow(length(results))
		  colors.graphlan = substr(colors.original, start = 1, stop = 7)
		  
		  for (result in results){
		    writestr = paste(result, "annotation",result,sep="\t")
		    write(writestr, 'annot.txt', append = TRUE)
		    milineage <- function(result, colors.graphlan) {
			    writestr = paste(result, "annotation_background_color", colors.graphlan[1],
			        sep = "\t")
		    }

		    write(writestr, 'annot.txt', append = TRUE)
		    writestr = paste(result, "clade_marker_shape\t*" ,sep="\t")
		    write(writestr, 'annot.txt', append = TRUE)
		    writestr = paste(result, "clade_marker_size\t20" ,sep="\t")
		    write(writestr, 'annot.txt', append = TRUE)
		    writestr = paste(result, "clade_marker_edge_color", colors.graphlan[1], sep="\t")
		    write(writestr, 'annot.txt', append = TRUE)
		    writestr = paste(result, "clade_marker_edge_width\t1" ,sep="\t")
		    write(writestr, 'annot.txt', append = TRUE)
		    write("", 'annot.txt', append = TRUE)

		    colors.graphlan = colors.graphlan[-1]
		  }
		  close(fileinfo)
		}
		## function for running graphlan

		run_graphlan = function(){
		  s1 = "export PATH=`pwd`/graphlan/:$PATH"
		  s2 = "export LC_ALL=en_US.UTF-8"
		  s3 = "export LANG=en_US.UTF-8"
		  s4 = "graphlan_annotate.py --annot annot.txt guide.txt guide.xml"
		  s5 = "graphlan.py guide.xml migraph.png --dpi 300 --size 3.5"
		  
		  call = paste(s1,s2,s3,s4,s5, sep = " && ")
		  system(call)
		}
		writetax(phyobj)
		writeannot(as.character(result$sig.lineage))
		run_graphlan()
	}

	# run
	result = milineage(phyobj, cov_n, cov_c, con_n, con_c, 
		n.resample = n.resample, fdr.alpha = fdr.alpha, min.depth = min.depth)
	# result$stacked_plot = get_stacked_plot(phyobj, result$sig.lineage)
	return(result)
	graphlan(phyobj, result)
	return(result)
}