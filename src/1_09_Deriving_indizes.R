#setting up working environment and load libraries
library(envimaR)
root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F:\\BEN\\edu")
source(file.path(root_folder, "mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))
#------------------------
#treepolygons
crown_seg <- readOGR(file.path(envrmt$path_clip, "aoi_crowns_clip.shp"))

crownH <- crown_seg$Height_m
crownA <- crown_seg$CA_m2
df <- data.frame(crownH = crownH, crownA = crownA)

boxplot(df$crownH, main = "Crown height")
boxplot(df$crownA, main = "Crown area")

#extern indizes
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Visible Atmospherically Resistant Index
VARI <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_VARI.tif"))
#test_crown <- readOGR(file.path(envrmt$path_lidar_processed, "crowns_1.shp"))

vari_ind <- extract(VARI, crown_seg, fun = mean)

#Triangular Greeness Index (to esatimate chlorophyll and (ind) plant nitrogen content)
TGI <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_TGI.tif"))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Area of AOI
area_total <- 0
for (i in length(shape1@polygons)) {
  area1<- shape1@polygons[[i]]@area
  area_total <- area_total + area1
}

#NA in CHM
chm_gaussian <- raster(file.path(envrmt$path_lidar_processed, "chm_aoi_gaussian.tif"))

values(chm_gaussian)
na_error <- (sum(is.na(values(chm_gaussian))) *100)/ length(values(chm_gaussian))
na_error

#
#l1 <- values(dsm_new)
#ldtm <- values(dtm_rst)
#l2 <- c()
#for (i in length(l1)) {
#  if(is.na(l1[i])){
#    l2[i] <- (ldtm[i])
#  } else{
#    l2[i] <- l1[i]
#  }
#}
#
#if(is.na(dsm_new[i])){
#  setValues(dsm_new[i], dtm_rst[i])
#}
