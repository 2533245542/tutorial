# This is a script for reproducing miLineage result

load("data_milineage_filtered.Rdata")
dt = phyobj

#*cova_var: selected categorical variable for 'covariates'
#*conf_var: selected categorical variable for 'confoundings'
#*cova_col: selected numerical columns of 'covariates'
#*conf_col: selected numerical columns of 'confoundings'

cova_var = NULL
cova_col = c("covariate")
conf_var = NULL
conf_col = c("confounding1","confounding2")

library(miLineage)
library(phyloseq)
source("functions_milineage_columninput.R")
source("functions_graphlan_annot.R")

# RUNNING MILINEAGE
results = run_QCAT(dt, cova_var, conf_var, cova_col, conf_col) 

# RETRIVE SIGNIFICANT LINEAGE
print(unlist(results$sig.lineage))

# WRITING GRAPHLAN INPUT FILES
writetax(dt) # taxonomy table (guide.txt)
writeannot(unlist(results$sig.lineage)) # annotation info (annot.txt)

# GENERATING VISUALIZED RESULT
run_graphlan()

# produce a phylogenetic tree with sig lineages indicated (migraph.png)
