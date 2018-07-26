## function for write tax
writetax = function(phy = GlobalPatterns){
  tax = tax_table(phy)
  ## unique tax
  tax = unique(tax)
  fileinfo = file('guide.txt')
  cat(file = fileinfo)
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

## function for running graphlan
run_graphlan = function(){
  s1 = "export PATH=`pwd`/graphlan/:$PATH"
  s2 = "export LC_ALL=en_US.UTF-8"
  s3 = "export LANG=en_US.UTF-8"
  s4 = "graphlan.py guide.txt migraph.png --dpi 300 --size 3.5"
  call = paste(s1,s2,s3,s4, sep = " && ")
  system(call)
}