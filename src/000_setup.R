#Skript for setting up a universal working environment on different computers
#Setting up project specific subfolders
#installation of devtools using devtools::install_github("envima/envimaR")


project_folders = c("data/",
                    "data/aerial/org/", "data/aerial/processed/","data/aerial/norm",
                    "data/lidar/org/", "data/lidar/clip/","data/lidar/processed/",
                    "data/lidar/test/","data/lidar/preprocessing/", "data/lidar/org/corrected/",
                    "data/lidar/preprocessing/level1",
                    "data/lidar/preprocessing/level2",
                    "data/lidar/preprocessing/level3",
                    "data/lidar/preprocessing/level4",
                    "data/lidar/preprocessing/level5",
                    "data/lidar/preprocessing/raster_level",
                    "data/lidar/preprocessing/raster_level/stdmetr",
                    "data/lidar/preprocessing/raster_level/shannon",
                    "data/lidar/preprocessing/raster_level/returns",
                    "data/grass/",
                    "data/height_models",
                    "data/data_mof/",
                    "data/tmp/",
                    "data/valid/",
                    "data/plots",
                    "data/species_det",
                    "run/",
                    "log/",
                    "mpg-envinfosys-teams-2018-bjcm_rs_18/src/",
                    "mpg-envinfosys-teams-2018-bjcm_rs_18/doc/",
                    "data/indices",
                    "data/species_comp")


# Set libraries
libs = c("lidR", "link2GI", "mapview", "raster", "rgdal","rlas", "sp")
lapply(libs, require, character.only = TRUE)

#Automatically set root direcory, folder structure and load libraries
envrmt = createEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd",
                   folders = project_folders,
                   path_prefix = "path_", libs = libs,
                   alt_env_id = "COMPUTERNAME", alt_env_value = "PCRZP",
                   alt_env_root_folder = "F:\\BEN\\edu")
