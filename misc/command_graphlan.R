s1 = "export PATH=`pwd`/graphlan/:$PATH"
s2 = "export LC_ALL=en_US.UTF-8"
s3 = "export LANG=en_US.UTF-8"
s4 = "graphlan.py guide.txt migraph.png --dpi 300 --size 3.5"
call = paste(s1,s2,s3,s4, sep = " && ")
# call
system(call)
library(png)
img <- readPNG(system.file("img", "Rlogo.png", package="png"))
img = readPNG("migraph.png")
str(img)
grid::grid.raster(img)

run_graphlan = function(){
  s1 = "export PATH=`pwd`/graphlan/:$PATH"
  s2 = "export LC_ALL=en_US.UTF-8"
  s3 = "export LANG=en_US.UTF-8"
  s4 = "graphlan.py guide.txt migraph.png --dpi 300 --size 3.5"
  call = paste(s1,s2,s3,s4, sep = " && ")
  system(call)
}

