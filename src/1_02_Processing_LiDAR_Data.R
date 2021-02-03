#Set working directory and load required packages
library(envimaR)
root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F:\\BEN\\edu")
source(file.path(root_folder, "mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))

#get all LAS files in a folder
las_files <- list.files(envrmt$path_lidar_org, pattern = glob2rx("*.las"),
                        full.names = TRUE)
#write LAX for each LAS
for (f in las_files) {
  writelax(f)
}

#read one single file
lidar_file <- readLAS(file.path(envrmt$path_lidar_org, "U4745630.las"))
plot(lidar_file, bg ="white", color = "Z")
rm(lidar_file)

#create catalog and set projection and parallel processing information
lcat <- catalog(envrmt$path_lidar_org)
opt_chunk_size(lcat) <- 500
lcat@output_options$output_files <- path.expand(envrmt$path_clip)

set_lidr_threads(4) #defauls is 2

#plot catalog
plot(lcat)

#clip catalog to the AOI
aoi <- readOGR(file.path(envrmt$path_data_mof, "AOI_mofshape.shp"))
aoi_bb <- bbox(aoi)#bounding box of shapefile


lasclipRectangle(lcat, xleft = aoi_bb[1], ybottom = aoi_bb[2],
                 xright = aoi_bb[3], ytop = aoi_bb[4])

#plotting controle
lcat_spatial <- as.spatial(lcat)
mapview(lcat_spatial)

        