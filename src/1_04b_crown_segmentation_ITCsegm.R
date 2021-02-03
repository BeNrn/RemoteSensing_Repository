#set working environment
library(envimaR)
library(itcSegment)
#library(uavRst)
root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F://BEN//edu")
source(file.path(root_folder, "mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))
#------------------------
chm <- raster::raster(file.path(envrmt$path_test, "chm_gaussian.tif"))#erstellt in CHM_raster.R

plot(chm)


#ITCsegment
crowns <- itcSegment::itcIMG(imagery = chm,
                             TRESHSeed = 0.45,
                             TRESHCrown = 0.55,
                             searchWinSize = 3,
                             epsg = 25832,
                             th = 2,
                             DIST = 40,
                             ischm = TRUE)
#alternative
#crowns <- uavRst::chmseg_ITC(chm = chm,
#                             EPSG = 25832,
#                             minTreeAlt = 2,
#                             TRESHSeed =0.45,
#                             TRESHCrown = 0.55,
#                             maxCrownArea = 40) #minTreeAlt mit Test optimieren
plot(crowns, add = TRUE)

#adding treetops

treetops <- tree_detection(chm, lmf(ws = 7, hmin = 2))#using local maximum filter: lidR:lmf
#ws: diameter of moving window
#hmin: minimum tree height
treetops@proj4string <- CRS("+init=epsg:25832")
plot(treetops, add = T, col = "black", legend = FALSE)
