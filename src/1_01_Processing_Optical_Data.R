#Set working directory and load required packages

library(raster)
# library(envimaR)
# root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd",
#                               alt_env_id = "COMPUTERNAME",
#                               alt_env_value = "PCRZP",
#                               alt_env_root_folder = "F:\\BEN\\edu")
# source(file.path(root_folder, "mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))
workingDir <- "D:/BEN/edu/data/"

#search images
#list_files <- list.files(file.path(envrmt$path_aerial_org), pattern = ".tif")
list_files <- list.files(paste0(workingDir, "aerial/org/"), pattern = ".tif")
list_files

#load images
imagelist <- list()
for (i in list_files) {
  imagelist[[i]] <- brick(paste0(workingDir, "aerial/org/", i))
}

#rename files in list to assign every image to its position
for (j in seq(1:length(imagelist))) {
  names(imagelist)[j] <- paste0(j, "_", list_files[j])
}

#projection
projct <- c()
for (x in imagelist) {
  v <- crs(x)
  projct <- c(projct, v)
}
projct

## projection the same?
if (length(unique(projct))==1){
  print("all projections are the same")
}else{
  print("at least one projection is different")
}

#load shapefile
shape1 <- readOGR(paste0(workingDir, "data_mof/", "AOI_mofshape.shp"),
                  layer = ogrListLayers(paste0(workingDir, "data_mof/", "AOI_mofshape.shp")))

#check coordinate system of the shapefile
crs(shape1)
#assign crs of the tifs
crs(shape1) <- projct[[1]]

#crop images on the extent of MOF
#the first two images are outside the bounding box
cropped <- lapply(imagelist[3:length(imagelist)], crop, shape1)

#mosaic overlapping images (IMPORTANT when using this for the first time -> set overwrite = T)
imagelist3_4 <- mosaic(cropped$`3_476000_5630000.tif`, cropped$`4_476000_5630000_1.tif`, fun = min, overwrite = T, filename = paste0(workingDir, "aerial/processed/", "image3_4.tif"))

imagelist5_6 <- mosaic(cropped$`5_476000_5632000.tif`, cropped$`6_476000_5632000_1.tif`, fun = min, overwrite = T, filename = paste0(workingDir, "aerial/processed/", "image5_6.tif"))

#merging of the images (IMPORTANT when using this for the first time -> set overwrite = T)
aoiimg <- merge(imagelist3_4, imagelist5_6, cropped$`7_478000_5630000.tif`, cropped$`8_478000_5632000.tif`, overwrite= T, filename=paste0(workingDir, "aerial/processed/", "AOIimage.tif"))

#plot of the optical data
plotRGB(brick(paste0(workingDir, "aerial/processed/", "AOIimage.tif")))
plot(shape1, add=TRUE)
