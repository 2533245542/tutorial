mistudio_file_storefile = function(name, path){ 
    # store user input file for reproducibilty; also rename to file to the correct name
    file.copy(path, "data")
    file_default_name = unlist(strsplit(path, "/"))[length(unlist(strsplit(path, "/")))]
    file.rename(paste0("data", "/", file_default_name), paste0("data", "/", name))
}
# store_file(input$file1$name, input$file1$datapath, "data")