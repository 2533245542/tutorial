# TODO test this function
mistudio_milineagec_compositional_reproduce = function(
	phyobj_original,
    phyobj_filtered,
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
	hidex
){
	clean_up_folder = function(){
		# each clean_up_folder function, though have the same name, is different from file to file
		# for example, this one is specific for this file
		if(file.exists("reproduce/reproduce_milineagec_compositional/phyobj.RData")){
			file.remove("reproduce/reproduce_milineagec_compositional/phyobj.RData")
		}
		if(file.exists("reproduce/reproduce_milineagec_compositional/reproduce_mistudio_milineagec_compositional.R")){
			file.remove("reproduce/reproduce_milineagec_compositional/reproduce_mistudio_milineagec_compositional.R")
		}
		if(file.exists("reproduce/reproduce_milineagec_compositional/last_image_before_download")){
			unlink("reproduce/reproduce_milineagec_compositional/last_image_before_download/*") # using unlnk to support the * operator
		}
		return()
	}
	write_phyobj = function(phyobj){
		save(phyobj, file = "reproduce/reproduce_milineagec_compositional/phyobj.RData")
		return()
	}
	write_reproduce_R = function(
		phyobj,
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
		create_file_with_header = function(path){
			file.create(path, showWarnings = TRUE)
			write("# load and install necessary packages", path, append = TRUE)
			write("source(\"scripts/installpkgs.R\")", path, append = TRUE)
			write("# load microbiome data", path, append = TRUE)
			write("load(\"phyobj.RData\")", path, append = TRUE)
			write("# load necesary functions", path, append = TRUE)
			write("source(\"scripts/functions.R\")", path, append = TRUE) 
			return(path)
		}
		write_parameter = function(
			path,
			phyobj,
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
			write_equation_parameter = function(var){
				name = substitute(var)
				content = var
			  	if(is.null(content)){
			   	 return(paste(name, "=", "NULL"))
			  	}			  
			  	if(is.numeric(content) || is.logical(content)){
			  		returnstring = paste(content, collapse = ",")
			  		returnstring = paste0("c(", returnstring, ")")
			  		returnstring = paste(name, "=", returnstring)
			  	} else {
			  		returnstring = paste(content, collapse = '","')
			  		returnstring = paste0("c(\"", returnstring, "\")")
			  		returnstring = paste(name, "=", returnstring)
			  	}
			  	
			  	return(returnstring)
			}			
			write("# paremeters", path, append = TRUE)
			write("phyobj = phyobj", path, append = TRUE)
			write(write_equation_parameter(arrange), path, append = TRUE)
			write(write_equation_parameter(lineage), path, append = TRUE)
			write(write_equation_parameter(stratify), path, append = TRUE)
			write(write_equation_parameter(substratify), path, append = TRUE)
			write(write_equation_parameter(continuous), path, append = TRUE)
			write(write_equation_parameter(categorical), path, append = TRUE)
			write(write_equation_parameter(subcategorical), path, append = TRUE)
			write(write_equation_parameter(xlab), path, append = TRUE)
			write(write_equation_parameter(ylab), path, append = TRUE)
			write(write_equation_parameter(palette), path, append = TRUE)
			write(write_equation_parameter(theme), path, append = TRUE)
			write(write_equation_parameter(hidex), path, append = TRUE)
			return(path)
		}
		write_get = function(path){
			write("# generate plot data", path, append = TRUE)
			write("data = get_plot_data(phyobj, arrange, lineage, stratify, substratify, continuous,categorical, subcategorical,xlab, ylab, palette, theme, hidex)", path, append = TRUE)
			return(path)
		}
		write_core_function = function(
			path,
			phyobj,
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
			write("# run function", path, append = TRUE)
			# make base plot
			if(arrange == "categorical"){
				if(is.null(xlab)){
					write("xlab = categorical", path, append = TRUE)
				}
				if(is.null(ylab)){
					write("ylab = \"relative abundance\"", path, append = TRUE)					
				}
				if(is.null(stratify)){ # no stratification
					if(is.null(categorical)){ # no sampling grouping
						write("p = ggplot(data, aes(x = sample_name, weight = rel, fill = name)) + geom_bar(position = \"fill\", width = 1, color = \"white\") + xlab(xlab) + ylab(ylab)", path, append = TRUE)	

					} else { # has sample grouping
						write("p = ggplot(data, aes(x = eval(as.symbol(categorical)), weight = rel_mean, fill = name)) + geom_bar(position = \"fill\", width = 1, color = \"white\") + xlab(xlab) + ylab(ylab)", path, append = TRUE)
					}
				}
				else if(length(stratify) == 1){ # one degree of stratification
					if(is.null(categorical)){ # no sample grouping
						write("p = ggplot(data, aes(x = sample_name, weight = rel, fill = name)) +	geom_bar(position = \"fill\", width = 1, color = \"white\") +facet_grid(rows = vars(eval(as.symbol(stratify)))) +xlab(xlab) + ylab(ylab)", path, append = TRUE)
						
					} else { # has sample grouping
						write("p = ggplot(data, aes(x = eval(as.symbol(categorical)), weight = rel_mean, fill = name)) +geom_bar(position = \"fill\", width = 1, color = \"white\") +facet_grid(rows = vars(eval(as.symbol(stratify)))) + xlab(xlab) + ylab(ylab)", path, append = TRUE)
						
					}
				} else if(length(stratify == 2)){ # two degrees of stratification
					if(is.null(categorical)){ # no sample grouping
						write("p = ggplot(data, aes(x = sample_name, weight = rel, fill = name)) +geom_bar(position = \"fill\", width = 1, color = \"white\") +facet_grid(rows = vars(eval(as.symbol(stratify[1]))), cols = vars(eval(as.symbol(stratify[2])))) + xlab(xlab) + ylab(ylab)", path, append = TRUE)
						
					} else { # has sample grouping
						write("p = ggplot(data, aes(x = eval(as.symbol(categorical)), weight = rel_mean, fill = name)) +geom_bar(position = \"fill\", width = 1, color = \"white\") +facet_grid(rows = vars(eval(as.symbol(stratify[1]))), cols = vars(eval(as.symbol(stratify[2])))) + xlab(xlab) + ylab(ylab)", path, append = TRUE)						
					}
				} else {
					warning("more than one stratify selected")
					return(NULL)
				}

			} else if(arrange == "continuous"){
				if(is.null(xlab)){
					write("xlab = continuous", path, append = TRUE)
				}
				if(is.null(ylab)){
					write("ylab = \"relative abundance\"", path, append = TRUE)
				}
				if(is.null(stratify)){
					write("p = ggplot(data, aes(x = eval(as.symbol(continuous)), weight = rel, fill = name)) +geom_histogram(binwidth = 1, position = \"fill\") + xlab(xlab) + ylab(ylab)", path, append = TRUE)
				}
				if(length(stratify) == 1){
					write("p = ggplot(data, aes(x = eval(as.symbol(continuous)), weight = rel, fill = name)) +geom_histogram(binwidth = 1, position = \"fill\") + facet_grid(rows = vars(eval(as.symbol(stratify[1])))) + xlab(xlab) + ylab(ylab)", path, append = TRUE)
					
				} else if(length(stratify) == 2){
					write("p = ggplot(data, aes(x = eval(as.symbol(continuous)), weight = rel, fill = name)) +geom_histogram(binwidth = 1, position = \"fill\") + facet_grid(rows = vars(eval(as.symbol(stratify[1]))), cols = vars(eval(as.symbol(stratify[2])))) + xlab(xlab) + ylab(ylab)", path, append = TRUE)				
				} else {
					warning("more than one stratify selected")
					return(NULL)
				}
			} else {
				warning("arrange is neither categorial nor continuous")
				return(NULL)
			}

			# add customized elements
			# steam line procedurally to add elements
			if(!is.null(palette)){			
				if(palette == "default"){
					p = p
				} else {
					if(palette == "1" || palette == "2" || palette == "3" || palette == "6" || palette == "7" || palette == "13"){
						write("palette = as.numeric(palette)", path, append = TRUE)
					}
					write("p = p + scale_fill_brewer(palette = palette)", path, append = TRUE)
				}			
			}                          
			if(!is.null(theme)){
				if(theme == "default"){
					p = p
				} else if(theme == "theme_bw"){
					write("p = p + theme_bw()", path, append = TRUE)
				} else if (theme == "theme_linedraw"){
					write("p = p + theme_light()", path, append = TRUE)
				} else if (theme == "theme_minimal"){
					write("p = p + theme_minimal()", path, append = TRUE)
				} else if (theme == "theme_classic"){
					write("p = p + theme_classic()", path, append = TRUE)
				} else if (theme == "theme_gray"){
					write("p = p + theme_gray()", path, append = TRUE)
				} else {
					warning("invalid milineaegc theme selected")
				}
			}
			if(hidex){
				write("p = p + theme(axis.text.x = element_blank(),	axis.ticks = element_blank())", path, append = TRUE)
				write("# if you want to adjust x ticks or y ticks more, check out this article:", path, append = TRUE)
				write("# http://www.sthda.com/english/wiki/ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels", path, append = TRUE)
			}
			# TODO: allow user to change legend as well
			write("p = p + labs(fill=\"lineage\")", path, append = TRUE)
			return(p)
		}

		"reproduce/reproduce_milineagec_compositional/reproduce_mistudio_milineagec_compositional.R" %>% 
		create_file_with_header(.) %>% 
		write_parameter(
			., 
			phyobj,
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
			hidex) %>% 
		write_get(.) %>% 
		write_core_function(
			., 
			phyobj,
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
			hidex)
		return()
	}
	# TODO finish write_last_image_before_download
	write_last_image_before_download = function(){
		return()
	}
	make_zip_file = function(){
		fs = dir(path = "reproduce/reproduce_milineagec_compositional/")
		fs = paste0("reproduce/reproduce_milineagec_compositional/", fs)
		fname = paste0("reproduce/reproduce_milineagec_compositional/", "reproduce_milineagec_compositional.zip")
		zip("reproduce/reproduce_milineagec_compositional/mistudio_milineagec_compositional_reproduce", fs)
	}
	clean_up_folder()
	write_phyobj(phyobj_filtered)
	write_reproduce_R(
	    phyobj_filtered,
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
		hidex
	)
	write_last_image_before_download()
	make_zip_file()
	return()
}