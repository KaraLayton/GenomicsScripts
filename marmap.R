##GEOGRAPHIC DISTANCE CALCULATION USING MARMAP
##first plot bathydata, and then use get.depth to select different locations and move points as necessary
library(marmap)
blues <- c("lightsteelblue4", "lightsteelblue3",
           "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))

##adjust lat and lons to zoom in on particular areas and move points
bathydata1 <- getNOAA.bathy(-54,-65,47,60,res=1,keep=T)
summary(bathydata1)
plot(bathydata1, image = TRUE, land = TRUE, lwd = 0.03,
     bpal = list(c(0, max(bathydata1), greys),
                 c(min(bathydata1), 0, blues)))
get.depth(bathydata1,distance=FALSE)
geodata2 <- read.csv("geosphere_dist_finaL_bioclim.csv", header=TRUE, row.names=1) #new geodata file created after having moved points

trans1 <- trans.mat(bathydata1,min.depth=10) ##can change min depth value
out_vst <- lc.dist(trans1, geodata2, res=c("dist"))
out_vst

##converting to matrix
out_vst <- as.matrix(out_vst)
write.csv(out_vst, "lcdist_bioclim.csv")
