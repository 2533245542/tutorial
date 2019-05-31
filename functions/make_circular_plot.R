# A function by https://rdrr.io/github/ying14/yingtools2/
# Converts taxonomy table to phylo object that handles singleton branches
myfun = function (x, data = parent.frame(), collapse.singles = FALSE, 
                  distinct.tree = TRUE, full.taxonomy.only = TRUE, ...) 
{
  err <- "Formula must be of the kind \"~A1/A2/.../An\"."
  if (length(x) != 2) 
    stop(err)
  if (x[[1]] != "~") 
    stop(err)
  f <- x[[2]]
  taxo <- list()
  while (length(f) == 3) {
    if (f[[1]] != "/") 
      stop(err)
    taxo[[deparse(f[[3]])]] <- data[[deparse(f[[3]])]]
    if (length(f) > 1) 
      f <- f[[2]]
  }
  taxo[[deparse(f)]] <- data[[deparse(f)]]
  taxo.data <- as.data.frame(taxo)
  if (distinct.tree) {
    taxo.data <- taxo.data %>% distinct()
  }
  if (full.taxonomy.only) {
    taxo.data <- taxo.data[!is.na(taxo.data[, 1]), ]
  }
  leaves.names <- as.character(taxo.data[, 1])
  taxo.data[, 1] <- 1:nrow(taxo.data)
  f.rec <- function(subtaxo) {
    u <- ncol(subtaxo)
    levels <- unique(subtaxo[, u])
    if (u == 1) {
      if (length(levels) != nrow(subtaxo)) 
        warning("Error, leaves names are not unique.")
      return(as.character(subtaxo[, 1]))
    }
    t <- character(length(levels))
    for (l in 1:length(levels)) {
      x <- f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 
                                                           1)])
      t[l] <- paste("(", paste(x, collapse = ","), ")", 
                    levels[l], sep = "")
    }
    return(t)
  }
  string <- paste("(", paste(f.rec(taxo.data), collapse = ","), 
                  ");", sep = "")
  phy <- read.newick(text = string)
  phy$edge.length <- rep(1, nrow(phy$edge))
  if (collapse.singles) {
    phy <- collapse.singles(phy)
  }
  phy$tip.label <- leaves.names[as.numeric(phy$tip.label)]
  return(phy)
}

library(ape)
library(phytools)
library(purrr)
data(carnivora)
t2 <- myfun(~SuperFamily/Family/Genus/Species, data=carnivora[,1:5])

par(lend=2)
plot(t2,edge.width=2,cex=0.6,no.margin=TRUE)
