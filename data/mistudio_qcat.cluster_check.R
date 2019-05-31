.Rarefy = function (otu.tab) 
{
  set.seed(100)
  depth = min(rowSums(otu.tab))
  otu.tab <- as.matrix(otu.tab)
  ind <- (rowSums(otu.tab) < depth)
  sam.discard <- rownames(otu.tab)[ind]
  otu.tab <- otu.tab[!ind, ]
  rarefy <- function(x, depth) {
    y <- sample(rep(1:length(x), x), depth)
    y.tab <- table(y)
    z <- numeric(length(x))
    z[as.numeric(names(y.tab))] <- y.tab
    z
  }
  otu.tab.rff <- t(apply(otu.tab, 1, rarefy, depth))
  rownames(otu.tab.rff) <- rownames(otu.tab)
  colnames(otu.tab.rff) <- colnames(otu.tab)
  return(otu.tab.rff)
}

library("car") # function recode
library("data.table")
library("miLineage")

load("~/Desktop/biome/_DATA/2018_1105_Sanjay/tax.st.Rdata")
load("~/Desktop/biome/_DATA/2018_1105_Sanjay/sanjay_st.RData")
load("~/Desktop/biome/_DATA/2018_1105_Sanjay/meta.st.Rdata")
load("~/Desktop/biome/_DATA/2018_1105_Sanjay/id.st.Rdata")
load("~/Desktop/biome/_DATA/2018_1105_Sanjay/otu.st.Rdata")

tax = tax.st[,-7]
otu = otu.st
id = id.st

# idx.otu = which(colSums(otu!=0)>=10 & substring(tax[,1], 4, ) =="Bacteria") #nothing is actually filtered
# otu = otu[,idx.otu]
# tax = tax[idx.otu,]

X1 = cbind(meta.st[,1])
X2 = cbind([,1], meta.st[,2])
X3 = cbind(meta.st[,1], meta.st[,2], meta.st[,1]*meta.st[,2])
X4 = cbind(meta.st[,1], meta.st[,2], meta.st[,2]*meta.st[,2])
X5 = cbind(meta.st[,1], meta.st[,2], meta.st[,1]*meta.st[,2],  meta.st[,2]*meta.st[,2], meta.st[,2]*meta.st[,2]*meta.st[,1] )

# adjust for baseline
base.idx = which(meta.st[,2]==1)
otu0 = otu[base.idx,]
otu2 = otu[-base.idx,]
set.seed(36)
otu2.rff = .Rarefy(otu2)
meta2.st = meta.st[-base.idx,]
meta2.st[,2] = meta2.st[,2]-2
id2 = id[-base.idx]
set.seed(36)
  otu0.rff = .Rarefy(otu0)
otu0.rff.prop = otu0.rff/rowSums(otu0.rff)
prop.base = NULL
for(i in 1:nrow(otu0)){
  prop.base = rbind(prop.base, rbind(otu0.rff.prop[i,], otu0.rff.prop[i,]) )
}


X1 = cbind(meta2.st[,1])
X2 = cbind(meta2.st[,1], meta2.st[,2])
X3 = cbind(meta2.st[,1], meta2.st[,2], meta2.st[,1]*meta2.st[,2])

set.seed(100)
otu2.rff = .Rarefy(otu2)
X = X2
res.qcat1.btw = QCAT.Cluster(id2, otu2, OTU.base=prop.base, X1, c(1), Tax = tax, n.perm  = 1)
# res.qcat1.all = QCAT.cluster(id2, "ALL", otu2, OTU.base=prop.base, X, c(1,3,5), Tax = tax, n.resample = 100)

sample_st = cbind(
  meta.st,
  id.st
)
rownames(sample_st) = rownames(otu)
a = sample_data(as.data.frame(sample_st, stringsAsFactors = FALSE))
b = otu_table(otu, taxa_are_rows = FALSE)
c = tax_table(tax)
d = phyloseq(a,b,c)
save(d, file = "st_timebase.RData")


