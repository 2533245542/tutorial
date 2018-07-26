# 
data("data.real")
tax = data.real$Tax.real
otu = data.real$OTU.real
sample = data.real$covariate.real
# otu.rff = .Rarefy(otu)
results.two = QCAT_GEE(otu.rff, sample, 1, sample, 1, tax, n.resample=100, fdr.alpha=0.05)
results.two$sig.lineage
# ...Converting to Phyloseq object
OTU = otu_table(otu, taxa_are_rows = FALSE)
TAX = tax_table(tax)
SAMPLE = sample_data(as.data.frame(sample))

physeq = phyloseq(OTU, TAX, SAMPLE)

# ...Conversion done
# ...Start extracting for QCAT_GEE
taxonomy = tax_table(physeq)
otutable = otu_table(physeq)
sampletable = sample_data(physeq)
dim(taxonomy)
dim(otutable)
dim(sampletable)

otu.rff = .Rarefy(otutable)
str(otu.rff)
str(sampletable)
str(taxonomy)
sampletable = as.matrix(sampletable)
taxonomy = as.matrix(taxonomy@.Data)
# ...Extraction done
# ...Start running QCAT_GEE
results.two = QCAT_GEE(otu.rff, sampletable, 1, sampletable, 1, taxonomy, n.resample=1000, fdr.alpha=0.05)
