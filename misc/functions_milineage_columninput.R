# Usage: results = run_QCAT(dt, covas, confs, covacols, confcols)

# Function for rarefy otu
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

##function for running QCAT
run_QCAT = function(dt = exampledata, covas, confs, covacols, confcols){
  library(caret)
  library(miLineage)

  # Get sample, otu, tax from physeq obj
  otu = otu_table(dt)
  tax = tax_table(dt)
  sample = sample_data(dt)

  # Handle dummy variable selection
  if(is.null(covas) & is.null(confs)){
    sample.QCAT_GEE = data.frame(matrix(ncol = 0, nrow = nrow(sample)))
    # No dummy conversion needed, do nothing
  } else  {

    # Yes dummy conversion needed,

    # Create full dummy variable matrix
    sample.data = sample@.Data
    sample.rownames = sample@row.names
    sample.colnames = sample@names
    sample.rebuild = as.data.frame(x = sample.data, row.names = sample.rownames, col.names = sample.colnames)
    sample.model = dummyVars(~., data = sample.rebuild)
    sample.dummy = predict(sample.model, sample.rebuild)

    # Select needed dummies
    conv_n_conf = c(covas, confs)
    indices.col = sapply(conv_n_conf, FUN = function(X)grep(X, colnames(sample.dummy)))
    sample.QCAT_GEE = sample.dummy[,indices.col]
  }

  # Handle numeric variable selection
  if(is.null(covacols) & is.null(confcols)){
    ifthisslot.isempty.itwill.printnull.whenentering.this = vector()
    # No numeric column selected, do nothing
  } else {

    # Yes numeric columns selected, cbind them to sample.QCAT_GEE
    sample.numcols = cbind(sample[,covacols], sample[,confcols]) #numeric variable matrix
    sample.QCAT_GEE = cbind(sample.QCAT_GEE, sample.numcols) #numeric variable + dummy variable matrix
  }

  set.seed(11) #for exampledata producibility

  otu.QCAT_GEE = .Rarefy(otu)
  tax.QCAT_GEE = tax@.Data
  colnames(tax.QCAT_GEE) = paste0("Rank", 1:ncol(tax.QCAT_GEE))
  sample.QCAT_GEE = as.matrix(sample.QCAT_GEE) # QCAT_GEE requires matrix form

  # Retrive cova's indices in sample.QCAT_GEE
  var.total = names(sample.QCAT_GEE)
  cova.total = c(covas, covacols)
  indices.cova.QCAT_GEE = sapply(cova.total, FUN = function(X) grep(X, cova.total))

  # otu Dimension Correction
  if(ncol(otu.QCAT_GEE) != nrow(tax.QCAT_GEE)){
    otu.QCAT_GEE = t(otu.QCAT_GEE)
  }


  QCAT_GEE(otu.QCAT_GEE, sample.QCAT_GEE, indices.cova.QCAT_GEE,
           sample.QCAT_GEE, indices.cova.QCAT_GEE, tax.QCAT_GEE,
           n.resample=1000, fdr.alpha=0.05)
}

# **********************testing1***************************
# data("GlobalPatterns")
# dt = GlobalPatterns
# covas = c("TGCGTT", "M31Fcsw", "Feces")
# confs = c("Soil", "SV1")
# covacols = NULL
# confcols = NULL
# **********************testing1***************************


# **********************testing2***************************
# data("exampledata")
# dt = exampledata
# covas = NULL
# confs = NULL
# covacols = "covariate"
# confcols = c("confounding1", "confounding2")
# **********************testing2***************************

# result = run_QCAT(dt,covas = covas,
#                   confs = confs, covacols = covacols,
#                   confcols = confcols)
