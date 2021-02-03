#LidR
#set working environment
library(envimaR)
library(uavRst)
root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F://BEN//edu")
source(file.path(root_folder, "mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))
#----------

#read CHM and tree position data

##uavRst##
chm <- raster::raster(file.path(envrmt$path_test, "chm_gaussian.tif"))#erstellt in CHM_with_gaussian.R
plot(chm)

#tree tops
treetops <- tree_detection(chm, lmf(ws = 7, hmin = 2))#using local maximum filter: lidR:lmf
#ws: diameter of moving window
#hmin: minimum tree height
treetops@proj4string <- CRS("+init=epsg:25832")
#plot(treetops, add = T, col = "black", pch = 1)

########################################
#Problem: treepos muss als Raster vorliegen, liegt aber als SpatialDataFrame vor.
########################################
# create spatial points data frame
spg <- data.frame(x = treetops@coords[,1], y = treetops@coords[2,], z = treetops$Z)
coordinates(spg) <- ~ x + y
# coerce to SpatialPixelsDataFrame
gridded(spg) <- TRUE
# coerce to raster
rasterDF <- raster(spg)
rasterDF

#see: python-kurs.github.io/
########################################



crownsRL <- chmseg_RL(treepos = rasterDF,
                      chm = chm,
                      maxCrownArea = 150,
                      exclusion = 0.2)
plot(crownsRL, add = T)
