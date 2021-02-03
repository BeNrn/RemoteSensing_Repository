#.libPaths("F:/lib")

#Set working directory and load required packages
require(uavRst)
require(lidR)
require(envimaR)
root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F:\\BEN\\edu")

source(file.path(root_folder, "/mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))
#-----------------------------------------------------
# define projection
proj4 = "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

######correct files have to be calculated once##########
#las_files = list.files(file.path(envrmt$path_lidar_org),
#                       pattern = glob2rx("*.las"), 
#                       full.names = TRUE)

#check extend of las files 
#uavRst::llas2llv0(las_files, paste0(envrmt$path_corrected, "/corrected"))
########################################################

#load corrected las files
las_files = list.files(envrmt$path_corrected, pattern = glob2rx("*.las"), 
                       full.names = TRUE)

#create catalog
mof_las <- uavRst::make_lidr_catalog(path = envrmt$path_corrected, 
                                     chunksize = 200, 
                                     chunkbuffer = 20, 
                                     proj4=proj4, 
                                     cores = 4)

#normalization
lidR::opt_output_files(mof_las) <- paste0(envrmt$path_preprocessing, sep = "/", "{ID}_normalized")

mof_norm <- lasnormalize(mof_las, tin())


#create new catalog with normalized files
norm_list <- list.files(envrmt$path_preprocessing, pattern = "normalized", full.names = TRUE)

mof_las_cor <- uavRst::make_lidr_catalog(norm_list, 
                                         chunksize = 200, 
                                         chunkbuffer = 20, 
                                         proj4=proj4, 
                                         cores = 4) 

#output options for ground classification
lidR::opt_output_files(mof_las_cor)<-paste0(envrmt$path_preprocessing, sep = "/", "{ID}_csf")

#define ground points
mof_las_cor<- lidR::lasground(mof_las_cor, csf())

las_list <- list.files(envrmt$path_preprocessing, pattern = "csf", full.names = TRUE)

#split LiDAR data in height layers
#for later call
#mof_las_cor <- uavRst::make_lidr_catalog(las_list, chunksize = 200, chunkbuffer = 20, proj4=proj4, cores = 4) 
#---------------

filteri = function(las,minZ = 0, maxZ = 5){
  las = readLAS(las)
  if (is.empty(lidR::lasfilter(las, Z >=minZ & Z < maxZ))) return(NULL)
  las = lidR::lasfilter(las, Z >=minZ & Z < maxZ)
  #grid = grid_metrics(las, res = 2, stdmetrics_i(Intensity))
  #grid = grid[[c(1,2,3,4)]] #
  return(las)
}

lidR::opt_output_files(mof_las_cor) <- paste0(envrmt$path_level1,"/{ID}_lv1_0_5m")
level1 = lidR::catalog_apply(mof_las_cor, filteri,0,5)
lidR:::catalog_laxindex(mof_las_cor)

lidR::opt_output_files(mof_las_cor) <- paste0(envrmt$path_level2,"/{ID}_lv2_5_10m")
level2 = lidR::catalog_apply(mof_las_cor, filteri,5,10)

lidR::opt_output_files(mof_las_cor) <- paste0(envrmt$path_level3,"/{ID}_lv3_10_15m")
level3 = lidR::catalog_apply(mof_las_cor, filteri,10,15)

lidR::opt_output_files(mof_las_cor) <- paste0(envrmt$path_level4,"/{ID}_lv4_15_20m")
level4 = lidR::catalog_apply(mof_las_cor, filteri,15,20)

lidR::opt_output_files(mof_las_cor) <- paste0(envrmt$path_level5,"/{ID}_lv5_20_50m")
level5 = lidR::catalog_apply(mof_las_cor, filteri,20,50)

#---------------------
#extract statistical metrics

#.stdmetrics_z
#[1] zmax
#[2] zmean
#[3] zsd
#[4] zskew
#[5] zkurt
#[6] zentropy
#[7] pzabovezmean
#[8] pzabove2
#[9 - 12] zq5 - 40

#standard metrics, see above
filterz = function(las){
  las = readLAS(las)
  if (is.empty(las)) return(NULL)
  grid = grid_metrics(las, res = 2, .stdmetrics_z)
  grid = grid[[c(1,2,3,6,7)]] 
  return(grid)
}

#shannon entropy
filt_sh = function(las){
  las = readLAS(las)
  if (is.empty(las)) return(NULL)
  grid_sh = grid_metrics(las, lidR::entropy(Z, by = 1), res = 2)
  return(grid_sh)
}

#returns
#filt_rt = function(las){
#  las = readLAS(las)
#  if (is.empty(las)) return(NULL)
#  grid_rt = grid_metrics(lidR::grid_density(las), res = 2)
#  return(grid_rt)
#}
#---------------------
#level 1 - stdmetrics
#---------------------
level1_list <- list.files(envrmt$path_level1, full.names = TRUE)

las_lv1_cat <- uavRst::make_lidr_catalog(level1_list, chunksize = 200, chunkbuffer = 20, proj4=proj4, cores = 4) 


lidR::opt_output_files(las_lv1_cat)<-paste0(envrmt$path_raster_level, sep = "/", "{ID}_zstats_lv1")
zstats = lidR::catalog_apply(las_lv1_cat, filterz)

#EXAMPLE for merge the generated tifs with included statistics
rastGes <- list.files(envrmt$path_raster_level, pattern = "zstats", full.names = TRUE)
for(i in seq(rastGes)){
  if(i == 1){
    rGes <- raster(rastGes[i])
  }else{
    rGes_tmp<- raster(rastGes[i])
    rGes <- merge(rGes, rGes_tmp)
  }
}
writeRaster(rGes, file.path(envrmt$path_tmp, "test_raster_lasfile.tif"))

#---------------------
#level 2 bis 5 - standard metrics
#---------------------
for(i in seq(2,5)){
  level_list <- list.files(file.path(paste0("~/edu/mpg-envinsys-plygrnd/data/lidar/preprocessing/level", i)), full.names = TRUE)
  las_lv_cat <- uavRst::make_lidr_catalog(level_list, chunksize = 200, chunkbuffer = 10, proj4=proj4, cores = 4)
  lidR::opt_output_files(las_lv_cat)<-paste0(envrmt$path_raster_level, "/{ID}_shannon_lv", i)
  shannon_lv = lidR::catalog_apply(las_lv_cat, filt_sh)
}


#merge raster images
for(j in seq(2,5)){
  rastGes <- list.files(file.path("F:/BEN/edu/data/lidar/preprocessing/raster_level/"), patter = paste0("shannon_lv", j) ,full.names = TRUE)
  for(i in seq(rastGes)){
    if(i == 1){
      rGes <- raster(rastGes[i])
    }else{
      rGes_tmp<- raster(rastGes[i])
      rGes <- merge(rGes, rGes_tmp)
    }
  }
  writeRaster(rGes, file.path(envrmt$path_lidar_processed, paste0("shannon_level", j, ".tif")))
}

#---------------------
#level 1 - shannon index
#---------------------
level1_list <- list.files(envrmt$path_level1, full.names = TRUE)

las_lv1_cat <- uavRst::make_lidr_catalog(level1_list, chunksize = 200, chunkbuffer = 10, proj4=proj4, cores = 4) 

lidR::opt_output_files(las_lv1_cat)<-paste0(envrmt$path_raster_level, sep = "/", "{ID}_shannon_lv1")
shannon_lv1 = lidR::catalog_apply(las_lv1_cat, filt_sh)
#merge raster images
rastGes <- list.files(envrmt$path_raster_level, pattern = "shannon_lv1", full.names = TRUE)
for(i in seq(rastGes)){
  if(i == 1){
    rGes <- raster(rastGes[i])
  }else{
    rGes_tmp<- raster(rastGes[i])
    rGes <- merge(rGes, rGes_tmp)
  }
}
writeRaster(rGes, file.path(envrmt$path_lidar_processed, "shannon_level1.tif"))

#---------------------
#level 2 bis 5 - shannon index
#---------------------
for(i in seq(2,5)){
  level_list <- list.files(file.path(paste0("~/edu/mpg-envinsys-plygrnd/data/lidar/preprocessing/level", i)), full.names = TRUE)
  las_lv_cat <- uavRst::make_lidr_catalog(level_list, chunksize = 200, chunkbuffer = 10, proj4=proj4, cores = 4)
  lidR::opt_output_files(las_lv_cat)<-paste0(envrmt$path_raster_level, "/{ID}_shannon_lv", i)
  shannon_lv = lidR::catalog_apply(las_lv_cat, filt_sh)
}


#merge raster images
for(j in seq(2,5)){
  rastGes <- list.files(file.path("F:/BEN/edu/data/lidar/preprocessing/raster_level/"), patter = paste0("shannon_lv", j) ,full.names = TRUE)
  for(i in seq(rastGes)){
    if(i == 1){
      rGes <- raster(rastGes[i])
    }else{
      rGes_tmp<- raster(rastGes[i])
      rGes <- merge(rGes, rGes_tmp)
    }
  }
  writeRaster(rGes, file.path(envrmt$path_lidar_processed, paste0("shannon_level", j, ".tif")))
}

#---------------------
#level 1 bis 5 - returns
#---------------------
#for(i in seq(1,5)){
#  level_list <- list.files(file.path(paste0("F:/BEN/edu/data/lidar/preprocessing/level", i)), full.names = TRUE)
#  las_lv_cat <- uavRst::make_lidr_catalog(level_list, chunksize = 200, chunkbuffer = 10, proj4=proj4, cores = 4)
#  lidR::opt_output_files(las_lv_cat) <- paste0(envrmt$path_returns, "/{ID}_returns_lv", i)
#  returns <- lidR::catalog_apply(las_lv_cat, filt_rt)
#}


#merge raster images
#for(j in seq(1,5)){
#  rastGes <- list.files(file.path("F:/BEN/edu/data/lidar/preprocessing/raster_level/returns"), patter = paste0("returns_lv", j) ,full.names = TRUE)
#  for(i in seq(rastGes)){
#    if(i == 1){
#      rGes <- raster(rastGes[i])
#    }else{
#      rGes_tmp<- raster(rastGes[i])
#      rGes <- merge(rGes, rGes_tmp)
#    }
#  }
#  writeRaster(rGes, file.path(envrmt$path_lidar_processed, paste0("returns_level", j, ".tif")))
#}
#
#--------------------------------------------------------

for(i in seq(level1_list)){
  if(i == 1){
    r_tmp <- grid_canopy(readLAS(level1_list[i]),
                         res = 2,
                         lidR::p2r(0.2,na.fill = knnidw()))
    r_tmp@crs <- CRS(proj4)
    rast_lv1 <- r_tmp
    print(paste0(i, ", 1"))
    
  }else{
    r_tmp <- grid_canopy(readLAS(level1_list[i]),
                         res = 2,
                         lidR::p2r(0.2,na.fill = knnidw()))
    r_tmp@crs <- CRS(proj4)
    rast_lv1 <- raster::merge(r_tmp, rast_lv1, overlap = TRUE) 
    print(paste0(i, ", 2"))
  }
}

writeRaster(rast_lv1, file.path(envrmt$path_raster_level, "rast_lv1.tif"), overwrite = TRUE)
plot(rast_lv1)
plot(grid_canopy(readLAS(file.path(envrmt$path_level1, "109_lv1_0_5m.las")),
                 res = 2,
                 lidR::p2r(0.2,na.fill = knnidw())),
     main = "Level 5")

rast_lv1 <- lapply(seq(1:10), function(i){
  if(i == 1){
    r_tmp <- grid_canopy(readLAS(level1[i]),
              res = 2,
              lidR::p2r(0.2,na.fill = knnidw()))
    r_tmp@crs <- CRS(proj4)
    rast_lv1 <- r_tmp
    class(rast_lv1)
    
  }else{
    r_tmp <- grid_canopy(readLAS(level1[i]),
                         res = 2,
                         lidR::p2r(0.2,na.fill = knnidw()))
    r_tmp@crs <- CRS(proj4)
    rast_lv1 <- raster::merge(r_tmp, rast_lv1) 
    print(paste(class(rast_lv1), "2"))
  }
  return(rast_lv1)
})

writeRaster(rast_lv1, file.path(envrmt$path_lidar_processed, "lev1_raster.tif"))

plot(grid_canopy(readLAS(file.path(envrmt$path_level5, "108_lv5_20_50m.las")),
                 res = 2,
                 lidR::p2r(0.2,na.fill = knnidw())),
     main = "Level 5")





lev1_list <- list.files(envrmt$path_level1, full.names = TRUE)

lev1_cat <- uavRst::make_lidr_catalog(lev1_list, 
                                      chunksize = 200, 
                                      chunkbuffer = 20, 
                                      proj4=proj4, 
                                      cores = 4) 
rasterfun = function(las){
  las = readLAS(las)
  if (is.empty(las)) return(NULL)
  grid <- grid_metrics(las, res = 2, .stdmetrics_i)
  return(grid)
}

lidR::opt_output_files(lev1_cat)<-paste0(envrmt$path_raster_level,"/{ID}_lv1rst_0_5m")
lev1_rst <- lidR::catalog_apply(lev1_cat, rasterfun)


for(i in 1:length(le)){
  las = readLAS(las_files[i])
  #las = lidR::lasfilter(las, Z >=0 & Z < 50)
  grid = grid_metrics(las, res = 2, .stdmetrics_i)
  crs(grid) =proj4                
  writeRaster(grid, filename = paste0("edu/mpg-envinsys-plygrnd/data/tmp/istats/",i,"_istats.tif"))
}





lidR::opt_output_files(mof_las_cor) <- paste0(envrmt$path_preprocessing, sep = "/", "{ID}_0_5")

filterz = function(las,minZ = 0, maxZ = 5){
  las = readLAS(las)
  if (is.empty(las)) return(NULL)
  las = lidR::lasfilter(las, Z >=minZ & Z < maxZ)
  grid = grid_metrics(las, res = 2, .stdmetrics_z)
  grid = grid[[c(1,2,3,6,7)]] #
  return(grid)
}




# statistics for total height
lidR::opt_output_files(mof_las_cor)<-paste0(envrmt$path_data_tmp,"zstats/{ID}_zstats")
zstats = lidR::catalog_apply(mof_snip_ground_csf, filterz,0,50)

lidR::opt_output_files(mof_snip_ground_csf)<-paste0(envrmt$path_data_tmp,"istats/{ID}_istats")
istats = lidR::catalog_apply(mof_snip_ground_csf, filteri,0,50)




lidR::opt_output_files(mof_snip_ground_csf)<-paste0(envrmt$path_data_tmp,"level3/{ID}_level3_istats")
level3 = lidR::catalog_apply(mof_snip_ground_csf, filteri,10,15)

lidR::opt_output_files(mof_snip_ground_csf)<-paste0(envrmt$path_data_tmp,"level4/{ID}_level4_istats")
level4 = lidR::catalog_apply(mof_snip_ground_csf, filteri,15,20)

lidR::opt_output_files(mof_snip_ground_csf)<-paste0(envrmt$path_data_tmp,"level5/{ID}_level5_istats")
level5 = lidR::catalog_apply(mof_snip_ground_csf, filteri,20,25)

lidR::opt_output_files(mof_snip_ground_csf)<-paste0(envrmt$path_data_tmp,"level6/{ID}_level6_istats")
level6 = lidR::catalog_apply(mof_snip_ground_csf, filteri,25,30)

lidR::opt_output_files(mof_snip_ground_csf)<-paste0(envrmt$path_data_tmp,"level7/{ID}_level7_istats")
level7 = lidR::catalog_apply(mof_snip_ground_csf, filteri,30,50)

#------------------------------------------------------------------------------------------------------
# create DSM
# add an output option FOR THE  pitfree algorithm
lidR::opt_output_files(mof_001)<-paste0(envrmt$path_height_models, sep = "/","{ID}_pfree_dsm")
dsm_pfree_csf <- lidR::grid_canopy(mof_001, 
                                   res = 2,
                                   lidR::pitfree(c(0,2,5,10,15), c(0, 0.5)))

# reclass spurious negative values
dsm_pfree_csf[dsm_pfree_csf<0]<-0

writeRaster(dsm_pfree_csf, filename = paste0(envrmt$path_data_lidar_prc,"dsm_pfree.tif"), overwrite = TRUE)

raster::plot(dsm_pitfree_csf,col=pal(32),main="csf pitfree c(0,2,5,10,15), c(0, 1.5)) 0.5 DSM")

#dtm using k means inverse distance approach
#k-nearest neighbour (KNN) and inverse-distance weighting (IDW) together
lidR::opt_output_files(mof_001_ground_csf)<-paste0(envrmt$path_species_det,"{ID}_knn_csf")
dtm_knn_csf = lidR::grid_terrain(mof_001_ground_csf, res=0.5,  algorithm = lidR::knnidw(k=50, p=3))
raster::plot(dtm_knn_csf,col=pal(32),main="csf knnidw terrain model")

# create chm

dsm_pfree_csf = resample(dsm_pfree_csf,dtm_knn_csf, method ="bilinear")
dsm_p2r_csf = resample(dsm_p2r_csf,dtm_knn_csf, method = "bilinear")
slope = terrain(dsm_pfree_csf,filename = paste0(envrmt$path_data_lidar_prc,"slope_pfree.tif"), opt = "slope", unit = "degrees", overwrite = TRUE)
writeRaster(dsm_pfree_csf, filename = paste0(envrmt$path_data_lidar_prc,"dsm_pfree.tif"), overwrite = TRUE)
writeRaster(dsm_p2r_csf, filename = paste0(envrmt$path_data_lidar_prc,"dsm_p2r.tif"), overwrite = TRUE)



#chm1 = dsm_pfree_csf - dtm_knn_csf
#chm2 = dsm_p2r_csf - dtm_knn_csf
#writeRaster(chm,filename = paste0(envrmt$path_data_lidar_prc,"pitfree-knn_chm.tif"), overwrite = TRUE)

#filter LAS catalog for different height levels

#### important: This only works if the following preparation steps have been conducted on LAS files:
#### 1. Find errors with uavRst::llas2llv0()
#### 2. Normalize Z values with lidR::lasnormalize()
#### 3. Reclassify ground returns with lidR::lasground()