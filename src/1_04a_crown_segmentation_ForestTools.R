#set working environment
library(envimaR)
library(ForestTools)
root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F://BEN//edu")
source(file.path(root_folder, "mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))
#------------------------
#load file !!CATALOG FEHLT!!
lidar_file <- readLAS(file.path(envrmt$path_test, "test_lasdata.las"))
lidar_file@crs <- CRS("+init=epsg:25832")

#create catalog and assign corrdinate system

#lcat <- catalog(envrmt$path_tmp)
#lcat@crs <- CRS("+init=epsg:25832")
#cores(lcat) <- 3L
#tiling_size(lcat) = 50
#plot(lcat)

#create dsm
dsm <-  grid_canopy(lidar_file, 0.5, p2r(subcircle = 0.2))
#dsm_rst <- as.raster(dsm)
dsm@crs <- CRS("+init=epsg:25832")

#create dtm
dtm <-  grid_terrain(lidar_file,
                     res = 0.5,
                     algorithm = kriging(k = 10))
#dtm_rst <- as.raster(dtm)
dtm@crs <- CRS("+init=epsg:25832")

#create chm
dsm_corr <- raster::resample(dsm, dtm , method = 'bilinear')
chm <- dsm_corr - dtm
plot(chm)

#ForestTools::vwf() treeposition
x <- function(x){x*0.05+0.6}#FUNCTION??
treetops <- vwf(CHM = chm, winFun = x, minHeight = 3)#minheight??

#alternative####################
# treetops <- tree_detection(chm, lmf(ws = 7, hmin = 2))#using local maximum filter: lidR:lmf
# #ws: diameter of moving window
# #hmin: minimum tree height
# treetops@proj4string <- CRS("+init=epsg:25832")
# plot(treetops, add = F)
################################
plot(treetops, add=TRUE, pch=10, col="red", cex=0.5)

#ForestTools::mcws() segmentation between treepos
crowns <- mcws(treetops = treetops, CHM = chm, format = "polygons", minHeight = 3, verbose = FALSE)
plot(crowns, add = TRUE)
