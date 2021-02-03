library(envimaR)
root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F:\\BEN\\edu")
source(file.path(root_folder, "mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))
source(file.path(envrmt$path_src, "Function_Validation.R"))
#-------------------------------------------------
#Die folgenden Änderungen benötigen folgende Anpassungen um das neue Skript zum laufen zu bringen:
#testval <- raster::raster(file.path(envrmt$path_tmp, "testalle.tif"),layer = ogrListLayers(paste(file.path(envrmt$path_tmp), "testalle.tif", sep = "/")))
#writeRaster(testval, file.path(envrmt$path_valid, "validierungsgebiet.tif"))
#-------------------------------------------------

#load testvalidaion (chm)
val_gebiet <- raster::raster(file.path(envrmt$path_valid, "validierungsgebiet.tif"))
val_gebiet@crs <- crs("+init=epsg:25832 +proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs ")

#filter
#chm_gaussian <- focal(testval, w = matrix(1, ncol = 3, nrow = 3), fun = mean)
chm_gaussian <- focal(val_gebiet, w = focalWeight(val_gebiet,d=c(1,1), type= "Gauss"))
chm_nofilter <- val_gebiet
#plot(chm_gaussian)

##load validation points
valid_pos <- readOGR(file.path(envrmt$path_valid, "Val_Tree_pos_Group.shp"),
                     layer = ogrListLayers(paste(envrmt$path_valid,  "Val_Tree_pos_Group.shp", sep = "/")))
###################################################################################################
#segmentation
#itcSegment
crowns <- uavRst::chmseg_ITC(chm = chm_gaussian,
                             EPSG = 25832,
                             movingWin = 3,
                             minTreeAlt = 12,
                             TRESHSeed =0.45,
                             TRESHCrown = 0.55,
                             maxCrownArea = 100)
#forest tools
x <- function(x){x*0.025+0.6}#FUNCTION??
treetops <- vwf(CHM = chm_gaussian, winFun = x, minHeight = 8)#minheight??

#ForestTools::mcws() segmentation between treepos
crowns <- mcws(treetops = treetops, CHM = chm_gaussian, format = "polygons", minHeight = 3, verbose = FALSE)

#lidR
treetops <- tree_detection(chm_gaussian, ws = 3, hmin = 8)#local maxima
treetops@crs <- CRS("+init=epsg:25832 +proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
#plot(treetops, add = T, col = "black", pch = 1)
crowns <- chmseg_RL(chm = chm_gaussian,
                      treepos = treetops,
                      maxCrownArea = 150,
                      exclusion = 0.2)
###################################################################################################

#transform coordinate system of validation points
vp <- spTransform(valid_pos,"+init=epsg:25832 +proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
identicalCRS(vp,crowns)

#assign a new id
vp$id <- as.numeric(vp$id)
for(i in seq(length(vp))){
  vp$id[i] <- i
}
#identify validation points that are assigned to more than one tree
#correction is the amount of correctly assigned points which share multiple crowns
correction <- c(0)

for(i in seq(length(vp@data$id))){
  singlePoint <- SpatialPointsDataFrame(coords = matrix(c(vp@coords[i,1], vp@coords[i,2]), nrow = 1, ncol = 2),
                                        proj4string = vp@proj4string,
                                        coords.nrs = vp@coords.nrs,
                                        data = data.frame(vp@data$id[i]))
  
  hit <- GISTools::poly.counts(singlePoint, crowns)
  if(sum(hit) > 1){
    correction <- correction + (sum(hit) - 1)
  }
}
print(correction)

#plot
plot(chm_gaussian, main = "ITC - Validierung")

plot(crowns, add = TRUE)

plot(vp, add = TRUE, pch = 1, col = "blue")

val_ft(vp = vp, x = crowns, corr = correction)
