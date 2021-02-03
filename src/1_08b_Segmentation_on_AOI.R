########latest changes
#buffer of spades split = 10, default = 0

#Set working directory and load required packages
#.libPaths("F:/lib")

library(envimaR)
library(itcSegment)
library(uavRst)
library(rLiDAR)
library(SpaDES.tools)
root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F:\\BEN\\edu")
source(file.path(root_folder, "/mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))
#---------------------------
#########################################
#assign coordinate system
proj4 = "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

#lidR alternative, lcat of AOI
lcat <- catalog(envrmt$path_lidar_org)
opt_cores(lcat) <- 3
opt_output_files(lcat) <- paste0(envrmt$path_clip, sep = "/", "las_aoi_clipped")
opt_chunk_size(lcat) <- 500
lcat@proj4string <- CRS(proj4)

#plot(lcat)

#clip catalog to the AOI bounding
#Output = one tile for the whole MOF
aoi <- readOGR(file.path(envrmt$path_data_mof, "AOI_mofshape.shp"))
aoi_bb <- bbox(aoi)#bounding box of shapefile
lcat <- lasclipRectangle(lcat, xleft = aoi_bb[1], ybottom = aoi_bb[2],
                         xright = aoi_bb[3], ytop = aoi_bb[4])

####################################################################################################
#add output options for the ground classification files we want to create 
lidR::opt_output_files(lcat)<-paste0(envrmt$path_clip,"ground_pts_aoi_csf")

#analyze ground points
ground_pts_aoi_csf<- lidR::lasground(lcat, lidR::csf())

#add an output option for pitfree algorithm
lidR::opt_output_files(ground_pts_aoi_csf)<-paste0(envrmt$path_height_models,"dsm_pitfree_csf")
dsm_pitfree_csf <- lidR::grid_canopy(ground_pts_aoi_csf, res = 0.5, 
                                     lidR::pitfree(c(0,2,5,10,15), c(0, 0.5)))

#reclass negative values
dsm_pitfree_csf[dsm_pitfree_csf < minValue(dsm_pitfree_csf)]<-minValue(dsm_pitfree_csf)

summary(dsm_pitfree_csf)

#plot dsm
raster::plot(dsm_pitfree_csf,main="csf pitfree c(0,2,5,10,15), c(0, 1.5)) 0.5 DSM")


#dsm
dsm <-  grid_canopy(lcat, res = 0.5, subcircle = 0.2)
dsm <-  grid_canopy(lcat, res = 0.5, subcircle = 0.2, na.fill = "kriging", k = 10L, model = NULL)
dsm_rst <- as.raster(dsm)
dsm_rst@crs <- CRS("+init=epsg:25832")
writeRaster(dsm_rst, filename = file.path(envrmt$path_lidar_processed, "dsm_aoi.tif"))

#dtm
dtm <-  grid_terrain(lcat, res = 0.5, method = "knnidw", k = 10L)
dtm_rst <- as.raster(dtm)
dtm_rst@crs <- CRS("+init=epsg:25832")
writeRaster(dtm_rst, filename = file.path(envrmt$path_lidar_processed, "dtm_aoi.tif"))

#chm
dsm_rst_corr <- raster::resample(dsm_rst, dtm_rst , method = 'bilinear')

#for later call:
#dsm_rst_corr <- raster(file.path(envrmt$path_lidar_processed, "dsm_aoi.tif"))
#dtm_rst <- raster(file.path(envrmt$path_lidar_processed, "dtm_aoi.tif"))

chm <- dsm_rst_corr - dtm_rst
writeRaster(chm, file.path(envrmt$path_lidar_processed, "chm_aoi.tif"), overwrite = F)

#gaussian filter, defined by the focalWeight function
chm_gaussian <- focal(chm, w = focalWeight(chm,d=c(1,1), type= "Gauss"), na.rm = TRUE)

#plot(chm_gaussian)
writeRaster(chm_gaussian, file.path(envrmt$path_lidar_processed, "chm_aoi_gaussian.tif"),overwrite = TRUE)
#rm(chm, dsm, dsm_rst, dsm_rst_corr, dtm, dtm_rst, aoi_bb, aoi) clear RAM
####################################################################################################

#Segmentation
#1. retile the chm
chm_gaussian #<- raster(file.path(envrmt$path_lidar_processed, "chm_aoi_gaussian.tif"))
chm_list <- list()
#64 matrix
chm_list <- SpaDES.tools::splitRaster(chm_gaussian, nx = 8, ny = 8, buffer = 20)
chm_list <- chm_list[10:64]
#2.start segmentation
lapply(seq(chm_list), function(l){
  crownseg <- uavRst::chmseg_ITC(chm = chm_list[[l]],
                               EPSG = 25832,
                               movingWin = 3,
                               minTreeAlt = 12,
                               TRESHSeed =0.45,
                               TRESHCrown = 0.55,
                               maxCrownArea = 80)
  
  m <- as.numeric(l)
  m <- m+9
  writeOGR(crownseg, paste0(envrmt$path_lidar_processed, "crowns_", m, ".shp"),
           driver = "ESRI Shapefile",
           layer = paste0("crowns_", m))
  print(Sys.time())
})
#list all proceeded shapes
crown_names <- list.files(envrmt$path_lidar_processed, pattern = ".shp")

#read all proceeded shapes
crown_list <- list()
for(i in crown_names){
  crown_list[[i]] <- readOGR(file.path(envrmt$path_lidar_processed, i))
}

#bind all shapefiles together
for (l in 1:length(crown_list)){
  tem <- crown_list[[l]]
  if (l==1){
  segm_aoi <- tem
  }
  else{
    segm_aoi <- rbind(segm_aoi, tem)
  }
}

#plot(segm_aoi)
writeOGR(segm_aoi, paste0("F:/BEN/edu/data/lidar/processed/", "aoi_crowns.shp"),
         driver = "ESRI Shapefile",
         layer = paste0("aoi_crowns.shp"))

plot(chm_list[[1]])
plot(crown_list[[1]], add = TRUE)
####################################################################################################