mistudio_milineage_reproduce = function(path,
	phyobj_name,
	data_path,

	filter_rank_selection,
	filter_rank,
	filter_samvars_selection,
	filter_samvars,
	filter_otu_aggregate,
	filter_otu_zerofraction,
	filter_otu_minsum,

	milineage_function,
	milineage_cova,
	milineage_conf,
	milineage_cova_extra,
	milineage_conf_extra,
	milineage_mindepth,
	milineage_nresample,
	milineage_fdralpha,
	milineage_ZI.LB,
	milineage_testtype,

	arrange, 
	lineage, 
	stratify, 
	substratify, 
	continuous,
	categorical, 
	subcategorical,
	xlab, 
	ylab, 
	palette, 
	theme, 
	hidex){
	get_file_name = function(phyobj_name, datapath){
		# search in all .biom or .BIOM files
		result1 = grep(".biom", list.files("data"), value = TRUE, ignore.case = TRUE)
		result2 = grep(phyobj_name, list.files("data"), value = TRUE, ignore.case = FALSE)
		file_name = intersect(result1, result2)
		# search find all .RData files
		file_name = lapply(grep(".RData", list.files("data"), value = TRUE, ignore.case = FALSE),
				function(x){
					# for each .RData file, find if phyobj_name is inside
					attach(paste0("data","/", x))
					exist = (length(grep(phyobj_name, ls(pos = 2))) != 0)
					detach()
					if(exist){
						return(x)
					} else {
						return(NULL)
					}
				}
		)
		file_name[sapply(file_name, is.null)] <- NULL
		file_name = unlist(file_name)

		if(length(file_name) == 0){
			warning("milineage reproduce fail: could not find the phyobj being used")
		}
		if(length(file_name) >= 2){
			warning("milineage reproduce: same phyloseq object names used in different .RData. Using the first one.")
			file_name = file_name[1]
		}
		file.copy(paste0("data","/", file_name), paste0(datapath, "/", file_name))
		return(file_name)
	}
	write_equation_parameter = function(path, var){
		name = substitute(var)
		content = var
	  	if(is.null(content)){
	   	 	returnstring = paste(name, "=", "NULL")
	  	} else if(is.numeric(content) || is.logical(content)){
	  		returnstring = paste(content, collapse = ",")
	  		returnstring = paste0("c(", returnstring, ")")
	  		returnstring = paste(name, "=", returnstring)
	  	} else {
	  		returnstring = paste(content, collapse = '","')
	  		returnstring = paste0("c(\"", returnstring, "\")")
	  		returnstring = paste(name, "=", returnstring)
	  	}
	  	write(returnstring, path, append = TRUE)
	}	
	# get_file_name(input$phyloseq_select, "reproduce/reproduce_milineage/data")
	create_file = function(path){
		if (file.exists(path)) {
			file.remove(path)
		}
		file.create(path)
		# also add functions
		file.copy("data/", "reproduce/milineage/", recursive = TRUE)
		file.copy("functions/", "reproduce/milineage/", recursive = TRUE)
	}
	write_header = function(path, main, subtitle = NULL){	
		write("################################################################################", path, append = TRUE)
		write("#", path, append = TRUE)
		write(paste0("# ", main), path, append = TRUE)
		if(!is.null(subtitle)){
			write(paste0("# ", subtitle), path, append = TRUE)
		}
		write("#", path, append = TRUE)
		write("################################################################################", path, append = TRUE)
	}
	write_intall_package = function(path){
		write("source(\"scripts/installpkgs.R\")",path,  append = TRUE)
		write("source(\"scripts/get_phyobj.R\")",path,  append = TRUE)
		write("source(\"scripts/filter_phyobj.R\")",path,  append = TRUE)
		write("source(\"functions/mistudio_milineage.R\")", path, append = TRUE)
		write("source(\"functions/mistudio_milineage_compositional.R\")", path, append = TRUE)
	}
	write_read_file = function(path, phyobj_name, data_path, file_name){
		write_equation_parameter(path, phyobj_name)
		write_equation_parameter(path, data_path)
		write_equation_parameter(path, file_name)
		write("phyobj = get_phyobj(phyobj_name, data_path, file_name)", path, append = TRUE)
		write("phyobj", path, append = TRUE)
	}
	write_filter_phyobj = function(path, 
		filter_rank_selection,
		filter_rank,
		filter_samvars_selection,
		filter_samvars,
		filter_otu_aggregate,
		filter_otu_zerofraction,
		filter_otu_minsum){
		browser()
		write_equation_parameter(path, filter_rank_selection)
		write_equation_parameter(path, filter_rank)
		write_equation_parameter(path, filter_samvars_selection)
		write_equation_parameter(path, filter_samvars)
		write_equation_parameter(path, filter_otu_aggregate)
		write_equation_parameter(path, filter_otu_zerofraction)
		write_equation_parameter(path, filter_otu_minsum)
		write("phyobj = filter_phyobj(\n    phyobj,\n    filter_rank_selection, \n    filter_rank, \n    filter_samvars_selection, \n    filter_samvars, \n    filter_otu_aggregate, \n    filter_otu_zerofraction, \n    filter_otu_minsum)", path, append = TRUE)
		write("phyobj", path, append = TRUE)
	}
	write_milineage = function(path,
		milineage_function,
	    milineage_cova,
	    milineage_conf,
	    milineage_cova_extra,
	    milineage_conf_extra,
	    milineage_mindepth,
	    milineage_nresample,
	    milineage_fdralpha,
	    milineage_ZI.LB,
	    milineage_testtype){
		write_equation_parameter(path, milineage_function)
		write_equation_parameter(path, milineage_cova)
		write_equation_parameter(path, milineage_conf)
		write_equation_parameter(path, milineage_cova_extra)
		write_equation_parameter(path, milineage_conf_extra)
		write_equation_parameter(path, milineage_mindepth)
		write_equation_parameter(path, milineage_nresample)
		write_equation_parameter(path, milineage_fdralpha)
		write_equation_parameter(path, milineage_ZI.LB)
		write_equation_parameter(path, milineage_testtype)
		write("milineage_result = mistudio_milineage(\n    phyobj, \n    milineage_function, \n    milineage_cova,\n    milineage_conf,\n    milineage_cova_extra,\n    milineage_conf_extra,\n    milineage_mindepth,\n    milineage_nresample,\n    milineage_fdralpha,\n    milineage_ZI.LB,\n    milineage_testtype)", path, append = TRUE)
		write("milineage_result", path, append = TRUE)
	}
	write_compositional = function(path, 
		arrange, 
		lineage, 
		stratify, 
		substratify, 
		continuous,
		categorical, 
		subcategorical,
		xlab, 
		ylab, 
		palette, 
		theme, 
		hidex){
		write_equation_parameter(path, arrange)
		write_equation_parameter(path, lineage)
		write_equation_parameter(path, stratify)
		write_equation_parameter(path, substratify)
		write_equation_parameter(path, continuous)
		write_equation_parameter(path, categorical)
		write_equation_parameter(path, subcategorical)
		write_equation_parameter(path, xlab)
		write_equation_parameter(path, ylab)
		write_equation_parameter(path, palette)
		write_equation_parameter(path, theme)
		write_equation_parameter(path, hidex)
		write("milineage_compositional = mistudio_milineage_compositional(\n    phyobj,\n    arrange, \n    lineage, \n    stratify, \n    substratify, \n    continuous,\n    categorical, \n    subcategorical,\n    xlab, \n    ylab, \n    palette, \n    theme, \n    hidex)", path, append = TRUE)
		write("milineage_compositional", path, append = TRUE)
	}
	make_zip_file = function(phyobj_name){
		zip_file_name = paste0(phyobj_name, "_milineage", "_reproduce")
		file_for_zipping = dir(path = paste0("reproduce/", "milineage"), full.names = TRUE)
		zip(zip_file_name, file_for_zipping)
	}
	
	file_name = get_file_name(phyobj_name, data_path)
	create_file(path)
	write_header(path, "Install/update missing packages and load functions")
	write_intall_package(path)
	write_header(path, "Read file and convert to phyloseq object")
	write_read_file(path, phyobj_name, data_path, file_name)
	write_header(path, "Filter phyloseq object")
	write_filter_phyobj(path, filter_rank_selection, filter_rank, filter_samvars_selection, filter_samvars, filter_otu_aggregate, filter_otu_zerofraction, filter_otu_minsum)
	if(milineage_function == "QCAT"){
		write_header(path, "Run miLineage", "Reference: A general framework for association analysis of microbial communities on a taxonomic tree. Tang ZZ, Chen G, Alekseyenko AV, Li H (2017). Bioinformatics, 33, 1278-1285.")
	} else if(milineage_function == "QCAT_GEE"){
		write_header(path, "Run miLineage", "Reference: A general framework for association analysis of microbial communities on a taxonomic tree. Tang ZZ, Chen G, Alekseyenko AV, Li H (2017). Bioinformatics, 33, 1278-1285.")
	} else if(milineage_function == "ZIGDM"){
		write_header(path, "Run miLineage", "Reference: Zero-inflated generalized Dirichlet multinomial regression model for microbiome compositional data analysis. Tang ZZ and Chen G (2018). Biostatistics, kxy025.")
	} else {
		warning("mistudio_milineage_reproduce: invalid milineage function")
	}
	write_milineage(path, milineage_function, milineage_cova,milineage_conf,milineage_cova_extra,milineage_conf_extra,milineage_mindepth,milineage_nresample,milineage_fdralpha,milineage_ZI.LB,milineage_testtype)
	write_header(path, "Make compositional plot")
	write_compositional(path, arrange, lineage, stratify, substratify, continuous,categorical, subcategorical,xlab, ylab, palette, theme, hidex)
	make_zip_file(phyobj_name)
}


# path = "reproduce/milineage/milineage_reproduce.R"
# phyobj_name = "milineagec_compositional_test"
# data_path = "data"
# file_name = "milineagec_compositional_test.RData"

# filter_rank_selection = "Phylum" 
# filter_rank = "Euryarchaeota"
# filter_samvars_selection = "SampleType"
# filter_samvars = c("Soil", "Skin", "Feces")
# filter_otu_aggregate = "Phylum"
# filter_otu_zerofraction = 0.9
# filter_otu_minsum = 20
# # TODO filter on mistudio web app is working; test reproduce to see if works
# milineage_function = "QCAT_GEE"
# milineage_cova = c("Sex", "Race")
# milineage_conf = NULL
# milineage_cova_extra = NULL 
# milineage_conf_extra = c("Color")
# milineage_mindepth = 0	
# milineage_nresample = 100
# milineage_fdralpha = 0.3
# milineage_ZI.LB = 10
# milineage_testtype = "PAR"

# arrange = "Categorical"
# lineage = "Ieradsrsad"
# stratify = "Sex"
# substratify = "Male"
# continuous = "Age"
# categorical = "Race"
# subcategorical = "White"
# xlab = "myxlab"
# ylab = "myylab"
# palette  = "Set1"
# theme = "blank"
# hidex = TRUE
# mistudio_milineage_reproduce(path,
# 	phyobj_name,
# 	data_path,

# 	filter_rank_selection,
# 	filter_rank,
# 	filter_samvars_selection,
# 	filter_samvars,
# 	filter_otu_aggregate,
# 	filter_otu_zerofraction,
# 	filter_otu_minsum,

# 	milineage_function,
# 	milineage_cova,
# 	milineage_conf,
# 	milineage_cova_extra,
# 	milineage_conf_extra,
# 	milineage_mindepth,
# 	milineage_nresample,
# 	milineage_fdralpha,
# 	milineage_ZI.LB,
# 	milineage_testtype,

# 	arrange, 
# 	lineage, 
# 	stratify, 
# 	substratify, 
# 	continuous,
# 	categorical, 
# 	subcategorical,
# 	xlab, 
# 	ylab, 
# 	palette, 
# 	theme, 
# 	hidex
# 	)