
# Prepare data
data(exampledata)
cova_var = c("21.61860")
cova_col = c("convariate")
conf_var = c("0.79515429")
conf_col = c("confounding1", "confounding2")

dt= exampledata

library(miLineage)
library(phyloseq)
source("functions_milineage_columninput.R")
source("functions_graphlan_annot.R")

# RUNNING MILINEAGE
results = run_QCAT(dt, cova_var, conf_var, cova_col, conf_col) 
          
          #*cova_var: selected categorical variable for 'covariates'
          #*conf_var: selected categorical variable for 'confoundings'
          #*cova_col: selected numerical columns of 'covariates'
          #*conf_col: selected numerical columns of 'confoundings'

# RETRIVE SIGNIFICANT LINEAGE
print(as.character(results$sig.lineage))

# WRITING GRAPHLAN INPUT FILES
writetax(phy = physeq()) # taxonomy table (guide.txt)
writeannot(as.character(results$sig.lineage)) # annotation info (annot.txt)

# GENERATING VISUALIZED RESULT
run_graphlan() 
          # produce a phylogenetic tree with sig lineages 
          # indicated (migraph.png)
