
# Usage: results = run_QCAT(GlobalPatterns, covas, confs)

#function edit rank from "Bacteria" to "k__Bacteria"
editrank = function(tax){
  for(rank in colnames(tax)){
    tax[,rank] = paste0(substr(tolower(rank), 1, 1), "__", tax[,rank])
  }
  colnames(tax) = paste0("Rank", seq(1:ncol(tax)))
  tax
}

#function for normalizing cova/conf
normalize = function(X1, samples){
  df = data.frame(c(1:nrow(X1)))
  for(sample in samples){
    #get col index
    col = which(X1 == sample, arr.ind = T)[1,2]
    #create normalize column
    df = cbind(df,ifelse(X1[, col] == sample, 1, 0) )
  }
  ##clean up and assign col names
  df = df[-1]
  colnames(df) = samples
  df
}

##function for count.rff
.Rarefy = function (otu.tab) 
{
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

## Append requested metadata columns to cova df or conf df
function(df, colnames){
  sample = sample_data(GlobalPatterns)
}
##function for running QCAT_GEE
run_QCAT = function(GlobalPatterns, covacols, confcols, covas, confs){
  # get X1, count1, tax1 from physeq()
  X1 = sample_data(GlobalPatterns)
  count1 = otu_table(GlobalPatterns)
  tax1 = tax_table(GlobalPatterns)
  
  # edit 
  X1 = X1[,-1]
  count1 = t(count1)
  tax1 = editrank(tax1)
  
  # normalize covas and confs
  df = normalize(X1, covas)
  df1 = normalize(X1, confs)
  X2 = cbind(df, df1)
  
  count.rff = .Rarefy(count1)
  
  # clean up 
  tax1 = tax1@.Data
  X2 = as.matrix(X2)
  # run
  QCAT_GEE(count.rff, X2, 1:length(df), X2, 1:length(df), tax1, n.resample=1, fdr.alpha=0.05 )
  
}








