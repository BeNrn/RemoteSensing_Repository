library(envimaR)
library(rgeos)
library(stringr)

root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F:\\BEN\\edu")
source(file.path(root_folder, "/mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))
#--------------------------------------------------------------------------------------
#######################################################################################
#neighbors (horizontal)
crowns <- readOGR(file.path(envrmt$path_clip, "PRobesample_Kronen.shp"))

AOI <- readOGR(file.path(envrmt$path_data_mof, "AOI_mofshape.shp"))

crowns_aoi <- raster::crop(crowns, AOI)
crowns_aoi$id <- 1:length(crowns_aoi)
#writeOGR(crowns_aoi, file.path(envrmt$path_species_comp, "crowns_aoi_ausschnitt.shp"), driver = "ESRI Shapefile", layer = "crowns_aoi_ausschnitt.shp")

centroids <- rgeos::gCentroid(crowns_aoi, byid = TRUE)

crowns_buffer <- rgeos::gBuffer(crowns_aoi, byid = TRUE, width = 15)

crowns_contain <- rgeos::gContains(crowns_buffer, centroids, byid = TRUE, returnDense = FALSE)

neighbor_p_tree <- data.frame(ID = crowns_aoi$id)

for(i in 1:8364){
  value <- as.numeric(length(crowns_contain[[i]]))
  neighbor_p_tree$neighbor[i] <- value 
}
saveRDS(neighbor_p_tree, file.path(envrmt$path_species_comp, "neighbors_15m.rds"))

neighbor_p_tree <- readRDS(file.path(envrmt$path_species_comp, "neighbors_15m.rds"))

#assign neighbors to polygons
crowns_aoi$Neighbors <- NA
for(s in seq(crowns_aoi$id)){
  crowns_aoi$Neighbors[s] <- neighbor_p_tree[neighbor_p_tree$ID == crowns_aoi$id[s],2]
  print(s)
}

writeOGR(crowns_aoi, file.path(envrmt$path_species_comp, "neighbors_15m.shp"), driver = "ESRI Shapefile", layer = "neighbors_15m.shp")
#######################################################################################
#vertical distribution (vertical)
raster_list <- list.files(file.path(envrmt$path_lidar_processed), pattern = "shannon", full.names = TRUE)
raster_list <- raster_list[which(stringr::str_detect(raster_list, "ovr") == FALSE)]
raster_list <- raster_list[which(stringr::str_detect(raster_list, "aux") == FALSE)]

allLvStack <- raster::stack(raster_list)

crowns_aoi <- readOGR(file.path(envrmt$path_species_comp, "crowns_aoi_ausschnitt.shp"))

#rasterite crown polys
ext = extent(crowns_aoi)
crown_rast = gdalUtils::gdal_rasterize(src_datasource = file.path(envrmt$path_species_comp, "crowns_aoi_ausschnitt.shp"), 
                                       dst_filename = file.path(envrmt$path_species_comp, "crowns_aoi_ausschnitt.tif"),
                                       l = "crowns_aoi_ausschnitt", tr = c(xres(allLvStack),yres(allLvStack)),
                                       te = c(ext[1],ext[3],ext[2],ext[4]),a = "id", output_Raster = T)

allLvStack <- raster::crop(allLvStack, ext)
allLvStack <- raster::resample(allLvStack, crown_rast)

allLvStack <- raster::stack(allLvStack[[1]], allLvStack[[2]], allLvStack[[3]], allLvStack[[4]], allLvStack[[5]])


#identify pixels with no tree polygon
vals = values(crown_rast)
index0 = which(vals == 0)

#assign id to polys
idPoly = vals[-index0]
df = data.frame(ID = idPoly)

#assign pixel values to df$id
for(i in 1:length(allLvStack@layers)){
  valsInd = values(allLvStack[[i]])
  valsInd = valsInd[-index0]
  name <- paste0("shn_", i)
  assign(name, valsInd)
}
rm(valsInd)

for(j in 1:5){
  name <- get(paste0("shn_", j))
  df = cbind(df, name)
}

names(df)[2:6] <- c("shannon_lv1","shannon_lv2", "shannon_lv3", "shannon_lv4", "shannon_lv5")

df[is.na(df)] <- 0

shannon <- data.frame(ID = NA, sh1 = NA, sh2 = NA, sh3 = NA, sh4 = NA, sh5 = NA)
for(k in 1:length(unique(df$ID))){
  sh_df <- df[df$ID == unique(df$ID)[k],]
  sh_lv1 <- mean(sh_df[[2]])
  sh_lv2 <- mean(sh_df[[3]])
  sh_lv3 <- mean(sh_df[[4]])
  sh_lv4 <- mean(sh_df[[5]])
  sh_lv5 <- mean(sh_df[[6]])
  shannon[k,] <- data.frame(ID = unique(df$ID)[k], sh1 = sh_lv1, sh2 = sh_lv2, sh3 = sh_lv3, sh4 = sh_lv4, sh5 = sh_lv5)
  print(k)
}


shannon$entropyPerMeter <- NA
for(d in 1:nrow(shannon)){
  if(shannon[d,6] != 0){
    shannon$entropyPerMeter[d] <- (((shannon[d,6]/30)*6) + shannon[d,5]/5 + shannon[d,4]/5+ shannon[d,3]/5 + shannon[d,2]/5)
  }else if(shannon[d,6] == 0 & shannon[d,5] != 0){
    shannon$entropyPerMeter[d] <- (shannon[d,5]/5 + shannon[d,4]/5+ shannon[d,3]/5 + shannon[d,2]/5)
  }else if(shannon[d,6] == 0 & shannon[d,5] == 0 & shannon[d,4] != 0){
    shannon$entropyPerMeter[d] <- (shannon[d,4]/5+ shannon[d,3]/5 + shannon[d,2]/5)
  }else if(shannon[d,6] == 0 & shannon[d,5] == 0 & shannon[d,4] == 0 & shannon[d,3] != 0){
    shannon$entropyPerMeter[d] <- (shannon[d,3]/5 + shannon[d,2]/5)
  }else if(shannon[d,6] == 0 & shannon[d,5] == 0 & shannon[d,4] == 0 & shannon[d,3] == 0 & shannon[d,2] != 0){
    shannon$entropyPerMeter[d] <- (shannon[d,2]/5)
  }else{
    shannon$entropyPerMeter[d] <- 0
  }
  print(d)
}


crowns_aoi$shannon <- NA

#assign species to polygons
for(s in seq(crowns_aoi$id)){
  if(is.na(unique(crowns_aoi$id[s] == shannon$ID)[2])){
    crowns_aoi$shannon[s] <- NA
  }else{
    crowns_aoi$shannon[s] <- shannon[shannon$ID == crowns_aoi$id[s],7]
  }
  print(s)
}

writeOGR(crowns_aoi, file.path(envrmt$path_species_comp, "shannonEntropy.shp"), driver = "ESRI Shapefile", layer = "shannonEntropy.shp")                    
#######################################################################################
#Höhe der Bäume (horizontal)
crowns_aoi <- readOGR(file.path(envrmt$path_species_comp, "crowns_aoi_ausschnitt.shp"))


