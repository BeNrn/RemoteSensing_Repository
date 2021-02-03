#Set working directory and load required packages
library(envimaR)
library(itcSegment)
library(uavRst)
root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F:\\BEN\\edu")
source(file.path(root_folder, "/mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))
#---------------------------
require(rLiDAR)
#load file
lidar_file <- readLAS(file.path(envrmt$path_clip , "las_mof_aoi_00001.las"))
writelax(lidar_file)

lidar_file@crs <- CRS("+init=epsg:25832")


#create catalog and assign corrdinate system
lcat <- catalog(file.path(envrmt$path_clip, "las_mof_aoi_00001.las"))
lcat@crs <- CRS("+init=epsg:25832")
cores(lcat) <- 3L
tiling_size(lcat) = 500#maximale Größe eines Catalog Fensters, wenn <=500, dann eine Kachel
#plot(lcat)

#create dsm
dsm <-  grid_canopy(lcat, res = 0.5, subcircle = 0.2)
dsm_rst <- as.raster(dsm)
dsm_rst@crs <- CRS("+init=epsg:25832")
#plot(dsm)

#create dtm
dtm <-  grid_terrain(lcat,
                     res = 0.5,
                     method = "kriging",
                     k = 10L)
dtm_rst <- as.raster(dtm)
dtm_rst@crs <- CRS("+init=epsg:25832")
#plot(dtm)

#create chm
dsm_rst_corr <- raster::resample(dsm_rst, dtm_rst , method = 'bilinear')
chm <- dsm_rst_corr - dtm_rst
writeRaster(chm, file.path(envrmt$path_tmp, "chm_testarea.tif"), overwrite = TRUE)
#plot(chm)

#require(graphics)
#plot(chm, col = grey.colors(12, start = 0, end = 0.9, gamma = 2.2, alpha = NULL))

# Load in AOI
#list.files(envrmt$path_tmp)
img_org <- raster::raster(file.path(envrmt$path_tmp, "chm_testarea.tif"))

#trinity <- matrix(1/9, ncol = 3, nrow = 3)

#gaussian filter, defined by the focalWeight function
chm_gaussian <- focal(img_org, w = focalWeight(img_org,d=c(1,1), type= "Gauss"))
#plot(chm_gaussian)
writeRaster(chm_gaussian, file.path(envrmt$path_tmp, "chm_gaussian.tif"),overwrite = TRUE)

#segmentation
crowns <- uavRst::chmseg_ITC(chm = chm_gaussian,
                             EPSG = 25832,
                             movingWin = 3,
                             minTreeAlt = 12,
                             TRESHSeed =0.45,
                             TRESHCrown = 0.55,
                             maxCrownArea = 80) 
plot(crowns, add = TRUE)
writeOGR(crowns, envrmt$path_log,
         driver = "ESRI Shapefile",
         layer = "crowns_aoi001")
