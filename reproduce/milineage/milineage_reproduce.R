################################################################################
#
# Install/update missing packages and load functions
#
################################################################################
source("scripts/installpkgs.R")
source("scripts/get_phyobj.R")
source("scripts/filter_phyobj.R")
source("functions/mistudio_milineage.R")
source("functions/mistudio_milineage_compositional.R")
################################################################################
#
# Read file and convert to phyloseq object
#
################################################################################
phyobj_name = c("milineagec_compositional_test")
data_path = c("data")
file_name = c("milineagec_compositional_test.RData")
phyobj = get_phyobj(phyobj_name, data_path, file_name)
phyobj
################################################################################
#
# Filter phyloseq object
#
################################################################################
filter_rank_selection = NULL
filter_rank = NULL
filter_samvars_selection = NULL
filter_samvars = NULL
filter_otu_aggregate = NULL
filter_otu_zerofraction = c(1)
filter_otu_minsum = c(0)
phyobj = filter_phyobj(
    phyobj,
    filter_rank_selection, 
    filter_rank, 
    filter_samvars_selection, 
    filter_samvars, 
    filter_otu_aggregate, 
    filter_otu_zerofraction, 
    filter_otu_minsum)
phyobj
################################################################################
#
# Run miLineage
# Reference: A general framework for association analysis of microbial communities on a taxonomic tree. Tang ZZ, Chen G, Alekseyenko AV, Li H (2017). Bioinformatics, 33, 1278-1285.
#
################################################################################
milineage_function = c("QCAT")
milineage_cova = c("Age%mistudio_seperator%Age","Sex%mistudio_seperator%Sex")
milineage_conf = c("Race%mistudio_seperator%Race","Color%mistudio_seperator%Color")
milineage_cova_extra = NULL
milineage_conf_extra = NULL
milineage_mindepth = c(0)
milineage_nresample = c(1000)
milineage_fdralpha = c(0.05)
milineage_ZI.LB = NULL
milineage_testtype = NULL
milineage_result = mistudio_milineage(
    phyobj, 
    milineage_function, 
    milineage_cova,
    milineage_conf,
    milineage_cova_extra,
    milineage_conf_extra,
    milineage_mindepth,
    milineage_nresample,
    milineage_fdralpha,
    milineage_ZI.LB,
    milineage_testtype)
milineage_result
################################################################################
#
# Make compositional plot
#
################################################################################
arrange = c("categorical")
lineage = c("Crenarchaeota")
stratify = c("Race","Color")
substratify = c("1%mistudio_seperator%Race","3%mistudio_seperator%Race","5%mistudio_seperator%Race","2%mistudio_seperator%Race","4%mistudio_seperator%Race","2%mistudio_seperator%Color","4%mistudio_seperator%Color","1%mistudio_seperator%Color","3%mistudio_seperator%Color")
continuous = NULL
categorical = c("Sex")
subcategorical = c("4%mistudio_seperator%Sex","3%mistudio_seperator%Sex","2%mistudio_seperator%Sex","1%mistudio_seperator%Sex")
xlab = c("szv")
ylab = c("pdf")
palette = c("Set2")
theme = c("theme_light")
hidex = c(TRUE)
milineage_compositional = mistudio_milineage_compositional(
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
milineage_compositional
