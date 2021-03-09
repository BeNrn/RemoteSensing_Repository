#script for tree crown segmentation using ForestTools 

#-------------------------------------------------------------------------------
#1 SET WORKING DIRECTORY AND LOAD PACKAGES
#-------------------------------------------------------------------------------
library(magrittr)
library(rLiDAR)
library(ForestTools)


#alternative using envimaR
# library(envimaR)
# root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
#                               alt_env_id = "COMPUTERNAME",
#                               alt_env_value = "PCRZP", 
#                               alt_env_root_folder = "F://BEN//edu")
# source(file.path(root_folder, "mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))

workingDir <- "D:/BEN/edu/data/"

#-------------------------------------------------------------------------------
#2 LOAD LIDAR FILE USING A CATALOG
#-------------------------------------------------------------------------------
#lidar_file <- paste0(workingDir, "lidar/test/test_lasdata.las") %>% readLAS()
#lidar_file@proj4string <- CRS("+init=epsg:25832")

#create catalog and assign corrdinate system
lcat <- paste0(workingDir, "lidar/test/test_lasdata.las") %>% catalog()
opt_chunk_size(lcat) <- 500
lcat@output_options$output_files <- path.expand(paste0(workingDir, "lidar/clip/"))
set_lidr_threads(4) #defauls is 2
opt_chunk_buffer(lcat) <- 10
#one tile, because it is smaller than 500m
plot(lcat)



#-------------------------------------------------------------------------------
#3 CREATE DSM USING A CATALOG
#-------------------------------------------------------------------------------
#set output location and name as 1dsm.tif
opt_output_files(lcat)<-paste0(workingDir, "lidar/test/{ID}dsm")
#create dsm  
dsm <-  grid_canopy(lcat, 
                    res = 0.5, 
                    algorithm = p2r(subcircle = 0.2))

#load dsm if already created
dsm <- raster(paste0(workingDir, "lidar/test/1dsm.tif"))

dsm@crs <- CRS("+init=epsg:25832")
plot(dsm)

#-------------------------------------------------------------------------------
#4 CREATE DTM USING A CATALOG
#-------------------------------------------------------------------------------
#set output location and name as 1dtm.tif
lidR::opt_output_files(lcat)<-paste0(workingDir, "lidar/test/{ID}dtm")
#create dtm
dtm <-  grid_terrain(lcat,
                     res = 0.5,
                     algorithm = kriging(k = 5))

#load dtm if already created
dtm <- raster(paste0(workingDir, "lidar/test/1dtm.tif"))

dtm@crs <- CRS("+init=epsg:25832")
plot(dtm)

#-------------------------------------------------------------------------------
#5 CREATE CHM
#-------------------------------------------------------------------------------
#create chm
#match the spatial resolution of the two tifs
dsm_corr <- raster::resample(dsm, dtm , method = 'bilinear')
chm <- dsm_corr - dtm
plot(chm)

#-------------------------------------------------------------------------------
#6 CROWN SEGMENTATION USING FOREST TOOLS
#-------------------------------------------------------------------------------
#problem with ForestTools: function has to be defined that determines the 
# variable moving window size. For example:
x <- function(x){x*0.05+0.6}
treetops <- ForestTools::vwf(CHM = chm, winFun = x, minHeight = 3)

plot(treetops, add=TRUE, pch=10, col="red", cex=0.5)

#ForestTools::mcws() segmentation between treepos
crowns <- mcws(treetops = treetops, CHM = chm, format = "polygons", minHeight = 3, verbose = FALSE)
plot(crowns, add = TRUE)
