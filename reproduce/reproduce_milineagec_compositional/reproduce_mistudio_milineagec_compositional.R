# load and install necessary packages
source("scripts/installpkgs.R")
# load microbiome data
load("phyobj.RData")
# load necesary functions
source("scripts/functions.R")
# paremeters
phyobj = phyobj
arrange = c("categorical")
lineage = NULL
stratify = NULL
substratify = NULL
continuous = NULL
categorical = NULL
subcategorical = NULL
xlab = NULL
ylab = NULL
palette = c("default")
theme = c("default")
hidex = c(FALSE)
# generate plot data
data = get_plot_data(phyobj, arrange, lineage, stratify, substratify, continuous,categorical, subcategorical,xlab, ylab, palette, theme, hidex)
# run function
xlab = categorical
ylab = "relative abundance"
p = ggplot(data, aes(x = sample_name, weight = rel, fill = name)) + geom_bar(position = "fill", width = 1, color = "white") + xlab(xlab) + ylab(ylab)
p = p + labs(fill="lineage")
