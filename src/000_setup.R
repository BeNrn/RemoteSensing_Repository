#--------NEUE VERSION
#installation über ---- devtools::install_github("envima/envimaR") ----
# Set project specific subfolders

project_folders = c("data/",                                 # data folders
                                        "data/aerial/org/", "data/aerial/processed/",
                                        "data/aerial/norm","data/lidar/org/", "data/lidar/clip/",
                                        "data/lidar/processed/", "data/lidar/test/","data/grass/", "data/height_models",
                                        "data/data_mof/", "data/tmp/", "data/valid/", "data/plots", "data/species_det",
                                        "data/lidar/preprocessing/",
                                        "run/", "log/", "data/lidar/org/corrected/",                         # bins and logging
                                        "mpg-envinfosys-teams-2018-bjcm_rs_18/src/",      # source code
                                        "mpg-envinfosys-teams-2018-bjcm_rs_18/doc/",
                                        "data/lidar/preprocessing/level1",
                                        "data/lidar/preprocessing/level2",
                                        "data/lidar/preprocessing/level3",
                                        "data/lidar/preprocessing/level4",
                                        "data/lidar/preprocessing/level5",
                                        "data/lidar/preprocessing/raster_level",
                                        "data/lidar/preprocessing/raster_level/stdmetr",
                                        "data/lidar/preprocessing/raster_level/shannon",
                                        "data/lidar/preprocessing/raster_level/returns",
                                        "data/indices", "data/species_comp")

#Check for needed packages
#required_packages <- c("envimaR", "lidR", "link2GI", "mapview", "raster", "rgdal","rlas", "sp",
#                       "RStoolbox", "rlist","graphics", "ForestTools", "itcSegment",
#                       "uavRst", "rLiDAR", "imager", "SpaDES.tools", "factoextra")


#if (length(setdiff(required_packages, rownames(installed.packages()))) > 0){
#  print(paste0(setdiff(required_packages, installed.packages()), "_package is not and will now be installed"))
#  install.packages(setdiff(required_packages, installed.packages()))
#} else {
#  print("all packages are installed")
#}

#rm(required_packages)

# Set libraries
libs = c("lidR", "link2GI", "mapview", "raster", "rgdal","rlas", "sp")
lapply(libs, require, character.only = TRUE)

#Automatically set root direcory, folder structure and load libraries
envrmt = createEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd",
                   folders = project_folders,
                   path_prefix = "path_", libs = libs,
                   alt_env_id = "COMPUTERNAME", alt_env_value = "PCRZP",
                   alt_env_root_folder = "F:\\BEN\\edu")
