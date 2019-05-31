 get_phyobj = function(phyobj_name, datapath, file_name){
 	last_char = tolower(tail(unlist(strsplit(file_name, "")), 1)) # get the last char of file_name
 	if(last_char == "a"){
 		env_temp = new.env()
 		load(paste0(datapath, "/", file_name), envir = env_temp)
 		phyobj = mget(phyobj_name, envir = env_temp)[[1]] #mget returns a list of objects but we are only using the first one and there is only one there	
 	} else if(last_char == "m"){
 		phyobj = import_biom(paste0(datapath, "/",file_name), silent = TRUE)
 	} else {
 		warning("get_phyobj: file_name is neither .Rdata or .BIOM")
 	}
 }
 # get_phyobj(input$phyloseq_select, "reproduce/reproduce_milineage/data", "asdf.BIOM")
 # filter_rank_selection
 # filter_rank
 # filter_samvars_selection
 # filter_samvars
 # filter_otu_aggregate
 # filter_otu_zerofraction
 # filter_otu_minsum


