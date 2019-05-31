mistudio_milineage_circular = function(phyobj, lineage){
	if(is.null(lineage)){
		return(NULL)
	}
	# A function by https://rdrr.io/github/ying14/yingtools2/
	# Converts taxonomy table to phylo object that handles singleton branches
	get_tax = function(phyobj){
		tax = tax_table(phyobj)@.Data
		return(tax)
	}

	myfun = function (x, data = parent.frame(), collapse.singles = FALSE, 
	                  distinct.tree = TRUE, full.taxonomy.only = TRUE, ...) 
	{
	  err <- "Formula must be of the kind \"~A1/A2/.../An\"."
	  if (length(x) != 2) 
	    stop(err)
	  if (x[[1]] != "~") 
	    stop(err)
	  f <- x[[2]]
	  taxo <- list()
	  while (length(f) == 3) {
	    if (f[[1]] != "/") 
	      stop(err)
	    taxo[[deparse(f[[3]])]] <- data[[deparse(f[[3]])]]
	    if (length(f) > 1) 
	      f <- f[[2]]
	  }
	  taxo[[deparse(f)]] <- data[[deparse(f)]]
	  taxo.data <- as.data.frame(taxo)
	  if (distinct.tree) {
	    taxo.data <- taxo.data %>% distinct()
	  }
	  if (full.taxonomy.only) {
	    taxo.data <- taxo.data[!is.na(taxo.data[, 1]), ]
	  }
	  leaves.names <- as.character(taxo.data[, 1])
	  taxo.data[, 1] <- 1:nrow(taxo.data)
	  f.rec <- function(subtaxo) {
	    u <- ncol(subtaxo)
	    levels <- unique(subtaxo[, u])
	    if (u == 1) {
	      if (length(levels) != nrow(subtaxo)) 
	        warning("Error, leaves names are not unique.")
	      return(as.character(subtaxo[, 1]))
	    }
	    t <- character(length(levels))
	    for (l in 1:length(levels)) {
	      x <- f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 
	                                                           1)])
	      t[l] <- paste("(", paste(x, collapse = ","), ")", 
	                    levels[l], sep = "")
	    }
	    return(t)
	  }
	  string <- paste("(", paste(f.rec(taxo.data), collapse = ","), 
	                  ");", sep = "")
	  phy <- read.newick(text = string)
	  phy$edge.length <- rep(1, nrow(phy$edge))
	  if (collapse.singles) {
	    phy <- collapse.singles(phy)
	  }
	  phy$tip.label <- leaves.names[as.numeric(phy$tip.label)]
	  return(phy)
	}

	library(ape)
	library(phytools)
	library(purrr)
	library(dplyr)
	library(ggtree)
	tax = get_tax(phyobj)
	tax = unique(tax)	
	for(i in 1:nrow(tax)){ # make sure there is no duplicate leaves(no NA)
	  for(j in 1:ncol(tax)){
	    if(is.na(tax[i,j])){
	      tax[i,j] = paste0(i,j)
	    }
	  }
	}
	t2 <- myfun(as.formula(paste0("~", paste(colnames(tax), collapse = "/"))), data= as.data.frame(tax))	
	p = ggtree(t2, layout="circular")
	color_index = sample(2:100, replace = FALSE)
	for(i in 1:length(lineage)){
	  tax_filtered = tax[as.logical(rowSums(tax == lineage[i], na.rm = TRUE)), ]
	  tax_filtered = tax_filtered[, ncol(tax_filtered)]
	  tax_filtered = as.character(tax_filtered) # obtain the tip node of the lineage chosen
	  tax_filtered = tax_filtered[!is.na(tax_filtered)] # remove NA if exists. NA in MRCA will not work.
	  node_num = MRCA(t2, tip=as.character(tax_filtered)) # put the tip node in MRCA to get node num of lineage
	  print(color_index[i])
    	p = p + geom_hilight(node=node_num, fill=color_index[i]) + 
        geom_cladelabel(node=node_num, label=lineage[i], align = TRUE, geom = "label") # use the node number to annotate    
	}
	return(p)
}








