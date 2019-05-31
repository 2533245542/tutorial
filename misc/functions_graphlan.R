function (OTU, Tree.file, gUniFrac.alpha = NULL, gUniFrac.rarefy = NULL, 
  pUniFrac.alpha = NULL, pUniFrac.rarefy = NULL, uwUniFrac.rarefy = NULL) 
{
  Tree = read.tree(Tree.file)
  if (!is.rooted(Tree)) {
    warning("The phylogenetic tree is not rooted. Rooting the tree using midpoint method and output to a new tree file.")
    Tree = midpoint(Tree)
    Tree.file = paste(Tree.file, "_rooted", sep = "")
    write.tree(Tree, file = Tree.file)
  }
  OTU = as.matrix(OTU)
  if (sum(!(colnames(OTU) %in% Tree$tip.label)) != 0) {
    warning("The OTU table contains OTUs that are not on the Tree. Removing those OTUs from the OTU table.")
    otu = OTU[, colnames(OTU) %in% Tree$tip.label]
  }
  else {
    otu = OTU
  }
  absent <- Tree$tip.label[!(Tree$tip.label %in% colnames(otu))]
  if (length(absent) != 0) {
    warning("The tree contains OTUs than are not in the OTU table. Trimming the tree to remove those OTUs and output to a new tree file.")
    Tree <- drop.tip(Tree, absent)
    Tree.file = paste(Tree.file, "_trimmed", sep = "")
    write.tree(Tree, file = Tree.file)
  }
  tmp = .RNewicktree(Tree.file)
  nbr = tmp$edgeCounts
  tree_tipLabel = tmp$tree_tipLabel
  tree_tipLabel_ID = tmp$tipId
  tree_upNode = tmp$edgeMatrix[, 1]
  tree_downNode = tmp$edgeMatrix[, 2]
  tree_brLen = tmp$edgeMatrix[, 3]
  otu = otu[, tree_tipLabel]
  n = nrow(otu)
  m = ncol(otu)
  if (is.null(rownames(otu))) {
    rownames(otu) <- paste("sample", 1:n, sep = "")
  }
  if (min(rowSums(otu)) <= 1) {
    otu_r = otu
  }
  else {
    set.seed(11)
    otu_r = .Rarefy(otu)
  }
  otu = otu/rowSums(otu)
  otu_r = otu_r/rowSums(otu_r)
  Dg = length(gUniFrac.alpha)
  Dg_type = gUniFrac.alpha
  Dg_norm = gUniFrac.rarefy
  if (Dg > 0 & length(Dg_norm) == 0) {
    Dg_norm = rep(1, Dg)
  }
  else if (Dg != length(Dg_norm)) {
    warning("Lengths of gUniFrac.alpha and gUniFrac.rarefy do not match. Do not rarefy data for generalized UniFrac calculations.")
    Dg_norm = rep(1, Dg)
  }
  Dp = length(pUniFrac.alpha)
  Dp_type = pUniFrac.alpha
  Dp_norm = pUniFrac.rarefy
  if (Dp > 0 & length(Dp_norm) == 0) {
    Dp_norm = rep(0, Dp)
  }
  else if (Dp != length(Dp_norm)) {
    warning("Lengths of pUniFrac.alpha and pUniFrac.rarefy do not match. Rarefy data for presence-weighted UniFrac calculations.")
    Dp_norm = rep(0, Dp)
  }
  Duw = length(uwUniFrac.rarefy)
  Duw_norm = uwUniFrac.rarefy
  res = .treeDist_v4(n, m, nbr, Dg, Dg_type, Dg_norm, Dp, 
    Dp_type, Dp_norm, Duw, Duw_norm, tree_tipLabel_ID, tree_brLen, 
    tree_upNode, tree_downNode, otu, otu_r)
  ndist = length(gUniFrac.alpha) + length(pUniFrac.alpha) + 
    length(uwUniFrac.rarefy)
  dist.name = NULL
  if (length(gUniFrac.alpha) > 0) {
    dist.name = c(dist.name, paste("gUniFrac", gUniFrac.alpha, 
      sep = "_"))
  }
  if (length(pUniFrac.alpha) > 0) {
    dist.name = c(dist.name, paste("pUniFrac", pUniFrac.alpha, 
      sep = "_"))
  }
  if (length(uwUniFrac.rarefy) > 0) {
    dist.name = c(dist.name, "uwUniFrac")
  }
  output <- array(NA, c(n, n, ndist), dimnames = list(rownames(otu), 
    rownames(otu), dist.name))
  start = 0
  end = 0
  index = 0
  if (length(gUniFrac.alpha) > 0) {
    for (k in 1:length(gUniFrac.alpha)) {
      index = index + 1
      start = end + 1
      end = end + n
      output[, , index] = res[start:end, ]
    }
  }
  if (length(pUniFrac.alpha) > 0) {
    for (k in 1:length(pUniFrac.alpha)) {
      index = index + 1
      start = end + 1
      end = end + n
      output[, , index] = res[start:end, ]
    }
  }
  if (length(uwUniFrac.rarefy) > 0) {
    index = index + 1
    start = end + 1
    end = end + n
    output[, , index] = res[start:end, ]
  }
  output.list = lapply(seq(dim(output)[3]), function(x) output[, 
    , x])
  names(output.list) = dimnames(output)[[3]]
  return(output.list)
}

Get_Dist = function(OTU, BC.rarefy = 0, Jaccard.rarefy = NULL) 
{
  otu <- as.matrix(OTU)
  n = nrow(otu)
  m = ncol(otu)
  if (is.null(rownames(otu))) {
    rownames(otu) <- paste("sample", 1:n, sep = "")
  }
  if (min(rowSums(otu)) <= 1) {
    otu_r = otu
  }
  else {
    set.seed(11)
    otu_r = .Rarefy(otu)
  }
  otu = otu/rowSums(otu)
  otu_r = otu_r/rowSums(otu_r)
  BC = length(BC.rarefy)
  Jaccard = length(Jaccard.rarefy)
  res = .Dist_v3(n, m, BC, BC.rarefy, Jaccard, Jaccard.rarefy, 
    otu, otu_r)
  # ndist = length(BC) + length(Jaccard)
  ndist = BC + Jaccard
  dist.name = NULL
  if (BC) {
    dist.name = c(dist.name, "Bray-Curtis")
  }

  if (Jaccard) {
    dist.name = c(dist.name, "Jaccard")
  }
  output <- array(NA, c(n, n, ndist), dimnames = list(rownames(otu), 
    rownames(otu), dist.name))
  start = 0
  end = 0
  index = 0
  if (BC) {
    index = index + 1
    start = end + 1
    end = end + n
    output[, , index] = res[start:end, ]
  }
  if (Jaccard) {
    index = index + 1
    start = end + 1
    end = end + n
    output[, , index] = res[start:end, ]
  }
  output.list = lapply(seq(dim(output)[3]), function(x) output[, 
    , x])
  names(output.list) = dimnames(output)[[3]]
  return(output.list)
}

## function for write tax
writetax = function(phy = GlobalPatterns){
  tax = tax_table(phy)
  ## unique tax
  tax = unique(tax)
  fileinfo = file('guide.txt')
  cat(file = fileinfo)
  i = 1
  while(i <= nrow(tax)){
    # parse NA
    names = as.character(tax[i,])
    index = as.logical(is.na(names))
    names = names[!index]
    ## set string
    writestr = paste0(names, collapse = ".")
    
    write(writestr, 'guide.txt', append = TRUE)
    i = i + 1
  }
  
  close(fileinfo)
  
}

## function for running graphlan
run_graphlan = function(){
  s1 = "export PATH=`pwd`/graphlan/:$PATH"
  s2 = "export LC_ALL=en_US.UTF-8"
  s3 = "export LANG=en_US.UTF-8"
  s4 = "graphlan.py guide.txt migraph.png --dpi 300 --size 3.5"
  call = paste(s1,s2,s3,s4, sep = " && ")
  system(call)
}

