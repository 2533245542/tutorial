provenence_milineage <- function(cova_var, cova_col, conf_var, conf_col, phyobj){
  # Turn string vector to string
  mipaste = function(varname, strings){
    if(is.null(strings)){
      return(paste(varname, "=", "NULL"))
    }
    
    returnstring = paste(strings, collapse = '","')
    returnstring = paste0("c(\"", returnstring, "\")")
    returnstring = paste(varname, "=", returnstring)
  }
  
  save(phyobj, file = "provenence_milineage/data_milineage_filtered.Rdata")
  
  fileinfo = file("provenence_milineage/script_miLineage.R")
  filename = "provenence_milineage/script_miLineage.R"
  cat(file = fileinfo)  
  
  write("# This is a script for reproducing miLineage result", filename, append = TRUE)
  write("", filename, append = TRUE)
  
  write("load(\"data_milineage_filtered.Rdata\")", filename, append = TRUE)
  write("dt = phyobj", filename, append = TRUE)
  write("", filename, append = TRUE)

  write("#*cova_var: selected categorical variable for 'covariates'", filename, append = TRUE)
  write("#*conf_var: selected categorical variable for 'confoundings'", filename, append = TRUE)
  write("#*cova_col: selected numerical columns of 'covariates'", filename, append = TRUE)
  write("#*conf_col: selected numerical columns of 'confoundings'", filename, append = TRUE)
  write("", filename, append = TRUE)
  
  write(mipaste("cova_var",cova_var), filename, append = TRUE)
  write(mipaste("cova_col",cova_col), filename, append = TRUE)
  write(mipaste("conf_var",conf_var), filename, append = TRUE)
  write(mipaste("conf_col",conf_col), filename, append = TRUE)
  write("", filename, append = TRUE)
  
  write("library(miLineage)", filename, append = TRUE)
  write("library(phyloseq)", filename, append = TRUE)
  write("source(\"functions_milineage_columninput.R\")", filename, append = TRUE)
  write("source(\"functions_graphlan_annot.R\")", filename, append = TRUE)
  write("", filename, append = TRUE)
  
  write("# RUNNING MILINEAGE", filename, append = TRUE)
  write("results = run_QCAT(dt, cova_var, conf_var, cova_col, conf_col) ", filename, append = TRUE)
  write("", filename, append = TRUE)
  
  write("# RETRIVE SIGNIFICANT LINEAGE", filename, append = TRUE)
  write("print(unlist(results$sig.lineage))", filename, append = TRUE)
  write("", filename, append = TRUE)
  
  write("# WRITING GRAPHLAN INPUT FILES", filename, append = TRUE)
  write("writetax(dt) # taxonomy table (guide.txt)", filename, append = TRUE)
  write("writeannot(unlist(results$sig.lineage)) # annotation info (annot.txt)", filename, append = TRUE)
  write("", filename, append = TRUE)
  
  write("# GENERATING VISUALIZED RESULT", filename, append = TRUE)
  write("run_graphlan()", filename, append = TRUE)
  write("", filename, append = TRUE)
  
  write("# produce a phylogenetic tree with sig lineages indicated (migraph.png)", filename, append = TRUE)
  close(fileinfo)
}

# # Prepare data
# cova_var = c("21.61860")
# cova_col = c("convariate")
# conf_var = c("0.79515429")
# conf_col = c("confounding1", "confounding2")
# provenence_milineage(cova_var, cova_col, conf_var, conf_col, exampledata)
