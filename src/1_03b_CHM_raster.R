#set working environment
library(envimaR)
root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F://BEN//edu")
source(file.path(root_folder, "mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))
#-------------------------------------------------------------------------------------

#load file
lidar_file <- readLAS(file.path(envrmt$path_test, "test_lasdata.las"))
lidar_file@proj4string <- CRS("+init=epsg:25832")

#create catalog and assign corrdinate system
lcat <- catalog(file.path(envrmt$path_test, "test_lasdata.las"))
opt_chunk_size(lcat) <- 500
lcat@output_options$output_files <- path.expand(envrmt$path_clip)

set_lidr_threads(4) #defauls is 2

opt_chunk_buffer(lcat) <- 10


#one tile, because it is smaller than 500m
plot(lcat)

#create dsm
dsm <-  grid_canopy(lidar_file, 0.5, p2r(subcircle = 0.2))
dsm@crs <- CRS("+init=epsg:25832")
plot(dsm)

#after first creation -> dtm is loaded below
#create dtm
#lidR::opt_output_files(lcat)<-paste0(envrmt$path_test,"/{ID}dtm")
#dtm <-  grid_terrain(lcat,
#                     res = 0.5,
#                     algorithm = kriging(k = 5))
dtm <- raster(paste0(envrmt$path_test,"/1dtm.tif"))

dtm@crs <- CRS("+init=epsg:25832")
plot(dtm)

#create chm
dsm_corr <- raster::resample(dsm, dtm , method = 'bilinear')
chm <- dsm_corr - dtm
#writeRaster(chm, file.path(envrmt$path_test, "chm_testarea.tif"))
plot(chm)

#-----------------------------------------------------------------------------

#gaussian filter, defined by the focalWeight function

chm_gaussian <- focal(chm, w = focalWeight(chm,d=c(1,1), type= "Gauss"))
plot(chm_gaussian)

writeRaster(chm_gaussian, file.path(envrmt$path_test, "chm_gaussian.tif"),overwrite = TRUE)
