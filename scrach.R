path = "gp_tax_uniquerow.txt"
write_tax = function(tax, path, removeNA = FALSE){
	if(file.exists(path)){
		file.remove(path)
	}
	for(i in 1:nrow(tax)){
		row_vector = tax[i, , drop = TRUE]
		row_string = paste0(row_vector, collapse = ".")
		if(removeNA){
			row_string = gsub(".NA", "", row_string)
		}
		write(row_string, path, append = TRUE)
	}
}             
