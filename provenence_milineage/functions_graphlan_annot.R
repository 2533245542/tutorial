
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
    writestr = paste(result, "annotation_background_color", colors.graphlan[1], sep="\t")
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
# data("exampledata")
# dt = exampledata
# covas = NULL
# confs = NULL
# covacols = "covariate"
# confcols = c("confounding1", "confounding2")
# # **********************testing2***************************
# 
# result = run_QCAT(exampledata,covas = covas,
#                   confs = confs, covacols = covacols,
#                   confcols = confcols)
# result = c("Erysipelotrichaceae", "Neisseriaceae", "Pseudomonadaceae")
# print("results done")
# print(as.character(result))
# writetax(dt)
# writeannot(as.character(result))
# print("writedone")
# run_graphlan()
# print("graphlan done")
# 
# file.remove("annot.txt", "guide.txt", "guide.xml", "migraph.png")
# file.remove("migraph.png")
