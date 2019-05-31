mistudio_miprofile_pca = function(
  phyobj= NULL, 
  gUniFrac.alpha = NULL, gUniFrac.rarefy = NULL, pUniFrac.alpha = NULL, pUniFrac.rarefy = NULL, uwUniFrac.rarefy = NULL,
  BC.rarefy = NULL, Jaccard.rarefy = NULL,
  tree.file = NULL,
  nperm = 1000,
  color = NULL,
  shape = NULL,
  facetrow = NULL,
  facetcol = NULL,
  label = NULL,
  label_size = NULL,
  label_x = NULL,
  label_y = NULL,
  label_rotate = NULL,
  palette = NULL,
  theme = NULL,
  point_size = NULL,
  point_opacity = NULL,
  xlab = NULL,
  ylab = NULL,
  width = NULL,
  height = NULL
  ){
  get_Dist.all = function(phyobj, 
    gUniFrac.alpha, gUniFrac.rarefy , pUniFrac.alpha , pUniFrac.rarefy, uwUniFrac.rarefy,
    BC.rarefy, Jaccard.rarefy,
    tree.file,
    nperm){
    stopifnot(!is.null(phyobj))
    adjust_otu_dimmension = function(otu, sample){
      if(dim(otu)[1] != dim(sample)[1]) {
        return(t(otu))
      } else {
        return(otu) 
      }
    }
    get_distance = function(){
      ensure_treefile = function(){
        # if user wants to calculate tree based distances, ensure tree file is available
        if(!is.null(gUniFrac.alpha) | !is.null(pUniFrac.alpha) | !is.null(uwUniFrac.rarefy)){
          if(!is.null(phy_tree(phyobj, errorIfNULL = FALSE))){
            write.tree(phy_tree(phyobj, "temp/miprofile_Get_Tree_Dist_treefile.tre"))
            Tree.file = "temp/miprofile_Get_Tree_Dist_treefile.tre" 
          } else{
            # if tree is not given in phyobj, it is better to have been uploaded by user now
            stopifnot(!is.null(tree.file))
          }     
        }     
      }

      if(!is.null(gUniFrac.alpha) || !is.null(pUniFrac.alpha) || !is.null(uwUniFrac.rarefy)){
        ensure_treefile()
        # tree based distances included
        stopifnot(!is.null(gUniFrac.rarefy) || !is.null(pUniFrac.rarefy))
        Dist.tree = Get_Tree_Dist(OTU = otu, Tree.file = tree.file, gUniFrac.alpha = gUniFrac.alpha, 
          gUniFrac.rarefy = gUniFrac.rarefy, pUniFrac.alpha = pUniFrac.alpha, 
          pUniFrac.rarefy = pUniFrac.rarefy, uwUniFrac.rarefy = uwUniFrac.rarefy ) 
        if(is.null(BC.rarefy) && is.null(Jaccard.rarefy)){
          Dist.all = Dist.tree
        } else {
          Dist.base = Get_Dist(OTU = otu, BC.rarefy = BC.rarefy, Jaccard.rarefy = Jaccard.rarefy) 
          Dist.all = c(Dist.base, Dist.tree)  
        }     
      } else{
        # tree based distances not included
        stopifnot(!is.null(BC.rarefy) || !is.null(Jaccard.rarefy))
        Dist.all = Get_Dist(OTU = otu, BC.rarefy = BC.rarefy, Jaccard.rarefy = Jaccard.rarefy)  
      }
      return(Dist.all)
    }
    # draw_all = function(dist.all = NULL, strata = NULL){
    #   stopifnot(!is.null(dist.all) || !is.null(strata))

    #   draw_one = function(dist.one, filename, strata){
    #     data = as.data.frame(cmdscale(as.dist(dist.one), k = 2))
    #     data = cbind(data, strata =strata)
    #     p = ggplot(data, aes(x = V1, y = V2, col = strata)) + 
    #     geom_point() + 
    #     labs(x  = "pc1", y = "pc2")   
    #     ggtitle(paste(filename, "PCA Plot")) +   
    #     theme(plot.title = element_text(hjust = 0.5))

    #     ggsave(filename = paste0("temp/",filename, ".pdf"), plot = p, device = 'pdf')
    #     return(p)
    #   }
    #   filenames = names(dist.all)

    #   graph = mapply(FUN = draw_one, dist.all, filenames, MoreArgs = list(strata  = strata))
    # }

    #sub-main
    otu = otu_table(phyobj)
    sample = sample_data(phyobj)
    otu = adjust_otu_dimmension(otu, sample)@.Data
    Dist.all = get_distance()
    return(Dist.all)
    # return(ggplot(mtcars, aes(x = mtcars$mpg, y = mtcars$cyl)) + geom_point())
  }

  draw_pca = function(
    Dist.all,
    phyobj,
    color,
    shape,
    facetrow,
    facetcol,
    label,
    label_size,
    label_x,
    label_y,
    label_rotate,
    palette,
    theme,
    point_size,
    point_opacity,
    xlab,
    ylab,
    width,
    height){
    get_sample = function(phyobj){
      sample = as.data.frame(sample_data(phyobj)@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
        col.names = sample_data(phyobj)@names, row.names = sample_data(phyobj)@row.names) # manually reassign column and row names    
      return(sample)
    }
    draw_one = function(
      Dist.one, 
      name,
      sample, 

      color,
      shape,
      facetrow,
      facetcol,
      label,
      label_size,
      label_x,
      label_y,
      label_rotate,
      palette,
      theme,
      point_size,
      point_opacity,
      xlab,
      ylab,
      width,
      height){
      data = as.data.frame(cmdscale(as.dist(Dist.one), k = 2))
      data = cbind(data, sample)

      p = ggplot(data, aes(x = V1, y = V2)) 
      # handle geom_point related elements: color, shape, point_size, point_opacity
      p = p + geom_point(aes_string(col = color, shape = shape), size = point_size , alpha = point_opacity)
      # handle geom_text related elements: label, label_size, label_x, label_y, label_rotate
      if(!is.null(label)){
        p = p + geom_text(aes_string(label = label), size = label_size, hjust = label_x, vjust = label_y, angle = label_rotate)        
      }
      
      # handle row/column faceting
      if(!is.null(facetrow) && !is.null(facetcol)){
        # handle row and column faceting
        p = p + facet_grid(cols = vars(eval(as.symbol(facetcol[1]))), rows = vars(eval(as.symbol(facetrow[1]))))  
      } else {
        if(!is.null(facetrow)){
          # handle row faceting
          p = p + facet_grid(rows = vars(eval(as.symbol(facetrow[1]))))  
        }        
        if(!is.null(facetcol)){
          # handle column faceting
          p = p + facet_grid(cols = vars(eval(as.symbol(facetcol[1]))))  
        }        
      }            
      # handle x, y lab customization
      if(xlab == ""){
        xlab = "PC1"
      }
      if(ylab == ""){
        ylab = "PC2"
      }
      p = p + xlab(xlab) + ylab(ylab)
      # handle palette and theme
      if(!is.null(palette)){      
        if(palette == "default"){
          p = p
        } else {
          if(palette == "1" || palette == "2" || palette == "3" || palette == "6" || palette == "7" || palette == "13"){
            palette = as.numeric(palette)
          }
          p = p + scale_fill_brewer(palette = palette)  
        }     
      }                          
      if(!is.null(theme)){
        if(theme == "theme_bw"){
          p = p
        } else if(theme == "theme_linedraw"){
          p = p + theme_linedraw()
        } else if (theme == "theme_light"){
          p = p + theme_light()
        } else if (theme == "theme_minimal"){
          p = p + theme_minimal()
        } else if (theme == "theme_classic"){
          p = p + theme_classic()
        } else if (theme == "theme_gray"){
          p = p + theme_gray()
        } else {
          warning("invalid milineaegc theme selected")
        }
      }
      # add the plot a title showing what distance it is
      p = p + ggtitle(paste(name, "PCA Plot")) + 
          labs(xlab = xlab, ylab = ylab) + 
          scale_alpha(guide = 'none') + scale_size(guide = 'none') 
          # + coord_fixed(ratio = 1)
      return(p)
      

      # p = ggplot(data, aes(x = V1, y = V2, col = strata)) + 
      #   geom_point() + 
      #   labs(x  = "pc1", y = "pc2") +
      #   ggtitle(paste(name, "PCA Plot")) +   # this adds the distance name
      #   theme(plot.title = element_text(hjust = 0.5))

      # # ggsave(filename = paste0("temp/",filename, ".pdf"), plot = p, device = 'pdf')
      # return(p)
    }
    sample = get_sample(phyobj)
    filenames = names(Dist.all)
    # browser()
    plot_list = mapply(FUN = draw_one, Dist.all, filenames, MoreArgs = list(
      sample,
      color,
      shape,
      facetrow,
      facetcol,
      label,
      label_size,
      label_x,
      label_y,
      label_rotate,
      palette,
      theme,
      point_size,
      point_opacity,
      xlab,
      ylab,
      width,
      height
      ), SIMPLIFY = FALSE
    )
    # plot = grid.arrange(plot_list[[1]],plot_list[[2]],ncol =1)
    plot = do.call(grid.arrange, c(plot_list, ncol = 1))
    return(plot)
  }

  # browser()
  plot = get_Dist.all(phyobj, 
    gUniFrac.alpha, gUniFrac.rarefy , pUniFrac.alpha , pUniFrac.rarefy, uwUniFrac.rarefy,
    BC.rarefy, Jaccard.rarefy,
    tree.file,
    nperm) %>% 
    draw_pca(
      ., 
      phyobj,
      color,
      shape,
      facetrow,
      facetcol,
      label,
      label_size,
      label_x,
      label_y,
      label_rotate,
      palette,
      theme,
      point_size,
      point_opacity,
      xlab,
      ylab,
      width,
      height)
    return(plot)

      
}

  

# get_pca_obj = function(mat = NULL, strata = NULL, name = NULL){
  
#   library(ade4)
  
#   stopifnot(!is.null(mat))
#   stopifnot(!is.null(strata))
#   stopifnot(!is.null(name))
  
#   strata = as.factor(strata)
#   quasi = quasieuclid(as.dist(mat))
#   pco = dudi.pco(quasi, scannf = F, nf = 2)
#   result = s.class(pco$li, fac = strata, 
#                    clabel=2, cpoint=1.5, grid=FALSE, addaxes=FALSE, 
#                    col=c( rgb(0,0,0,alpha=1),rgb(1,0,0,alpha=0.5) ), 
#                    sub=name, csub=3)
#   return(result)
# }

# pdf(file="figPC_status_normal.pdf", width=15, height=11, paper = "special")  
# par(mfrow = c(2,2))

# get_pca_obj(mat = mat, strata = strata, name = "Bray-Curtis") # name is the legend's name

# dev.off()

# function(){

# }