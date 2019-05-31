
get_pca_obj = function(mat = NULL, strata = NULL, name = NULL){
  
  library(ade4)
  
  stopifnot(!is.null(mat))
  stopifnot(!is.null(strata))
  stopifnot(!is.null(name))
  
  strata = as.factor(strata)
  quasi = quasieuclid(as.dist(mat))
  pco = dudi.pco(quasi, scannf = F, nf = 2)
  result = s.class(pco$li, fac = strata, 
                   clabel=2, cpoint=1.5, grid=FALSE, addaxes=FALSE, 
                   col=c( rgb(0,0,0,alpha=1),rgb(1,0,0,alpha=0.5) ), 
                   sub=name, csub=3)
  return(result)
}

pdf(file="figPC_status_normal.pdf", width=15, height=11, paper = "special")  
par(mfrow = c(2,2))

get_pca_obj(mat = mat, strata = strata, name = "Bray-Curtis")

dev.off()