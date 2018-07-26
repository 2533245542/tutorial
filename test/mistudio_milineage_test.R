source('../functions/mistudio_milineage.R')
library(phyloseq)

load('../data/exampledata.RData')
load("/Users/wzhou87/Desktop/data.toy_v3_test.RData")
data('GlobalPatterns')

# debug(mistudio_milineage)
mistudio_milineage(phyobj = data.toy_v3_test, cov_n = 'case', 
	n.resample = 1000, fdr.alpha = 0.05)
# mistudio_milineage(exampledata)


